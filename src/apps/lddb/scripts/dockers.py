
from abc import ABCMeta, abstractmethod, abstractproperty
import os
import random
import math
import cStringIO
import generator
import subprocess
import lddb_parsers
import json
import docking_pipeline
import shutil

class AbstractDocker(generator.AbstractGenerator):
  __metaclass__ = ABCMeta

  # Actual data generation
  @abstractmethod
  def run(self, input_params):
    raise NotImplementedError( "Should have implemented this" )

  # Copy data to an appropriate location on a file system and/or load data onto a database
  @abstractmethod
  def load(self, conn, d):
    raise NotImplementedError( "Should have implemented this" )

  # data locations to be synchronized to reduce data movement per job
  @abstractmethod
  def get_sync_info(job_dict):
    raise NotImplementedError( "Should have implemented this" )

  # expected output files.. the deployment will wait for this one
  @abstractmethod
  def get_output_fns(job_dict):
    raise NotImplementedError( "Should have implemented this" )

  # Copy data to a local running location and performs a light computation e.g., sed command
  @abstractmethod
  def preprocess(d, local_paths):
    raise NotImplementedError( "Should have implemented this" )

#=====================================================================================
class Dock6Docker(AbstractDocker):

  # In this case the run method calls a pipleline run method instead of performing
  # the actual computation itself. A pipeline is used when the computation
  # contains multiple steps and check pointing is required
  def run(self, d):
    pl = docking_pipeline.Dock6Pipeline(d)
    status = pl.run() or 'normal'
    output_prefix   = os.path.join(d['resource_prefix'],
                                   d['io_dir'],
                                   d['session_dir'],
                                   d['deployment_dir'],
                                   d['run_dir'])

    status_fn = str(os.path.join(output_prefix, 'status.txt'))
    with open(status_fn, 'w') as ofp:
      ofp.write(status)

  # Calling the pipeline preprocess method for consistency
  @staticmethod
  def preprocess(d, local_paths):
    docking_pipeline.Dock6Pipeline.preprocess(d, local_paths)

  @staticmethod
  def get_output_fns(d=None):
    out_dict = { 'status_fn': 'status.txt'
               }
    return out_dict

  # We might want to include receptor frame grid diretories here later on
  @staticmethod
  def get_sync_info(d):
    return []


  # Copy the grid output to the receptor frame directory if it's not there already
  # Copy the entire run directory to the job directory and load the scores onto the database
  def load(self, conn, prefix, d, local_paths):
    #trj_dir = os.path.join(d['dest_dir'], "recp{receptor_id}_model{model_id}/trj{trj_id}".format(**d))
    #frame_dir     = os.path.join(recpmodel_dir, "frame{t}".format(**d))
    frame_dir = d['frame_dir']

    grid_fn = 'grid.out'
    selected_spheres_fn = 'selected_spheres.sph'
    boxpdb_fn = 'step2_box.pdb'
    bmp_fn = 'step3_grid_score.bmp'
    nrg_fn = 'step3_grid_score.nrg'
    if not (os.path.exists(os.path.join(frame_dir, grid_fn)) and \
          os.path.exists(os.path.join(frame_dir, selected_spheres_fn)) and \
          os.path.exists(os.path.join(frame_dir, bmp_fn)) and \
          os.path.exists(os.path.join(frame_dir, nrg_fn)) and \
          os.path.exists(os.path.join(frame_dir, boxpdb_fn)) ):
      cmd_st = 'cp *grid* *.sph *box* {0}'.format(frame_dir)
      subprocess.call(cmd_st, shell = True, cwd = prefix, stdout=subprocess.PIPE,
                      stdin=subprocess.PIPE)
      print 'Copied all to %s' % frame_dir

    # Move the docked file to the unique name
    original_docked_fn = 'docked_ligand_scored.mol2'
    docked_fn = '{0}_level{1}_docked_{2}.mol2'.format(d['docker'], d['docker_level'], d['compoundstate_id'])
    cmd_st = 'mv {0} {1}'.format(original_docked_fn, docked_fn)
    subprocess.call(cmd_st, shell = True, cwd = prefix, stdout=subprocess.PIPE,
                      stdin=subprocess.PIPE)

    # Copy the docked file to the frame dir
    cmd_st = 'cp {0} {1}'.format(docked_fn, frame_dir)
    subprocess.call(cmd_st, shell = True, cwd = prefix, stdout=subprocess.PIPE,
                      stdin=subprocess.PIPE)
 
    d['dest_dir']        = os.path.join(local_paths['resource_prefix'],
                                        d['dbname'])
    job_dir = os.path.join(d['dest_dir'], 'job'+str(d['job_id']))
    if not os.path.exists(job_dir):
      os.makedirs(job_dir)

    cmd_st = 'cp * {0}'.format(job_dir)
    subprocess.call(cmd_st, shell = True, cwd = prefix, stdout=subprocess.PIPE,
                    stdin=subprocess.PIPE)
    #with open(os.path.join(job_dir, 'job_dict.txt'), 'w') as ofp:
    #  ofp.write(str(d))

    score_fn = os.path.join(job_dir, 'step4_dock6.out')

    results = lddb_parsers.parse_dock6(score_fn)
    
    st_out = cStringIO.StringIO()
    for r in results:
      r = [d['job_id']] + r
      r = map(str, r)
      st_out.write('\t'.join(r) + '\n')

    st_out.seek(0)

    cur = conn.cursor()
    cur.copy_from(st_out, 'Level1_dock6scores', columns = ('job_id', 'gs', 'gvdw', 'ges', 'inen'))
    cur.close()
    conn.commit()
    return bool(results), job_dir



#=====================================================================================

class FredDocker(AbstractDocker):
  #print '\t\t FredDocker made\n'

  #KTW@staticmethod
  #def runFred(job_dict):
  def run(self, job_dict):
    status = 'normal'###pl.run() or 'normal'
    output_prefix   = os.path.join(job_dict['resource_prefix'],
                                   job_dict['io_dir'],
                                   job_dict['session_dir'],
                                   job_dict['deployment_dir'],
                                   job_dict['run_dir'])
    status_fn = str(os.path.join(output_prefix, 'status.txt'))
    with open(status_fn, 'w') as ofp:
      ofp.write(status)

    job_dict['sdf_fn']       = os.path.splitext(job_dict['smi_fn'])[0] + '.sdf'

    cmd_st = ('nice -n 19   {omega2_bin} -in {smi_fn} -out {sdf_fn} ' +\
                           '-maxconfs 0 -maxtime 1200 -rms 0.25 -ewindow 20'
             ).format(**job_dict)
    #job_dict['recepmod_fn']      = "{resource_prefix}/{recepmod_fred_dir}/{recepmod_fred_fn}".format(**job_dict)
#
    run_at = "{output_prefix}".format(**job_dict)
    subprocess.Popen(cmd_st, cwd = run_at,
                     shell=True, stdout=subprocess.PIPE,
                     stdin=subprocess.PIPE).stdout.read()

    job_dict['recepmod_fn'] = os.path.splitext('protein.pdb')[0] + '.oeb.gz'
    if not (os.path.exists(os.path.join(output_prefix, job_dict['recepmod_fn']))):
      #prepare the receptor pdb into oeb.gz
      cmd_st = ('{apopdb2receptor} -pdb protein.pdb -site_residue PHE104 -solvent_residues NoSolvents -receptor {recepmod_fn}').format(**dict(job_dict.items()))
      subprocess.Popen(cmd_st, cwd = run_at,
                       shell=True, stdout=subprocess.PIPE,
                       stdin=subprocess.PIPE).stdout.read()


    #run fred docking
    job_dict['docked_fn'] = '{0}_level{1}_docked_{2}.sdf'.format(job_dict['docker'], job_dict['docker_level'], job_dict['compoundstate_id'])
    job_dict['score_fn']  = 'score.txt'
    job_dict['report_fn']  = 'report.txt'
    job_dict['settings_fn']  = 'settings.txt'
    cmd_st = ('nice -n 19 /damsl/projects/molecules/software/DD/OpenEye/centos/openeye/bin/fred -receptor {recepmod_fn} -dbase {sdf_fn} -dock_resolution High -docked_molecule_file {docked_fn} -score_file {score_fn} -report_file {report_fn} -settings_file {settings_fn} -save_component_scores true -hitlist_size 0').format(**dict(job_dict.items()))
    subprocess.Popen(cmd_st, cwd = run_at,
                       shell=True, stdout=subprocess.PIPE,
                       stdin=subprocess.PIPE).stdout.read()

  @staticmethod
  def preprocess(job_dict, local_paths):
    #print "\tFred preprocess"
    #print "\t\tJob_dict\n\t\t%s\n" % job_dict
    #print "\t\tLocal_paths\n\t\t%s\n" % local_paths
    #pasteing KTW

    job_dict['dest_dir']        = os.path.join(local_paths['resource_prefix'],
                                        job_dict['dbname'])
    frame_dir = os.path.join(job_dict['dest_dir'], "recp{receptor_id}/model{model_id}/trj{trj_id}/frame{t}".format(**job_dict))
    #frame_dir     = os.path.join(recpmodel_dir, "frame{t}".format(**d))
    job_dict['frame_dir'] = frame_dir

    #print 'frame_dir is: %s' % frame_dir
    if not os.path.exists(frame_dir):
      os.makedirs(frame_dir) 
      #print 'made %s ' % frame_dir
    frame_fn = os.path.join(frame_dir, 'protein.pdb')
    #copy pdb to the frame_dir
    if not os.path.exists(frame_fn):
      shutil.copy2(os.path.join(job_dict['trj_dir'], job_dict['pdb_fn']), frame_fn)
    '''
    #KTW - really need to define output_prefix?
    output_prefix   = os.path.join(job_dict['resource_prefix'],
                                   job_dict['io_dir'],
                                   job_dict['session_dir'],
                                   job_dict['deployment_dir'],
                                   job_dict['run_dir'])
    job_dict['output_prefix'] = output_prefix
    
    '''
    output_prefix = job_dict['output_prefix']
    #create the run_dir
    if not os.path.exists(output_prefix):
      os.makedirs(output_prefix)

    #copy only pdb from frame_dir to run_dir
    frame_fn = os.path.join(frame_dir, 'protein.pdb')
    cmd = 'cp {0} {1}'.format(frame_fn, os.path.join(output_prefix, 'protein.pdb'))
    subprocess.call(cmd, shell=True)

    #if it exists, copy oeb.gz from frame_dir to run_dir
    recmod_fn = os.path.join(frame_dir, 'protein.oeb.gz')
    cmd = 'cp {0} {1}'.format(recmod_fn, output_prefix)
    subprocess.call(cmd, shell=True)
   
    #copy license from frame_dir to run_dir
    cmd = 'cp {0} {1}'.format(local_paths['OE_LICENSE'], output_prefix)
    subprocess.call(cmd, shell=True)

    #write the smile file
    smile_fn = 'ligand.smi'
    job_dict['smi_fn'] = smile_fn
    with open(os.path.join(output_prefix, smile_fn), "w") as smi:
      smi.write(job_dict['smiles'])

  @staticmethod
  def get_output_fns(job_dict):
    #basename  = 'c{c_id}_m{m_id}_r{r_id}_t{taut}'.format(**job_dict)
    out_dict = { 'score_fn':  'score.txt'
                ,'docked_fn': '{0}_level{1}_docked_{2}.sdf'.format(job_dict['docker'], job_dict['docker_level'], job_dict['compoundstate_id'])
                ,'report_fn': 'report.txt'
                ,'status_fn': 'status.txt'
               }
    return out_dict
#
  @staticmethod
  def get_sync_info(job_dict):
    return [(job_dict['recepmod_fred_dir'], "*.oeb.gz")]

  # Copy oeb receptor if it isn't there, docked_fn, and score_fn to frame_dir
  def load(self, conn, prefix, job_dict, local_paths):
    
    frame_dir = job_dict['frame_dir']
    job_dict['recepmod_fn'] = 'protein.oeb.gz'

    if not (os.path.exists(os.path.join(frame_dir, job_dict['recepmod_fn']))):
      cmd_st = 'cp {0} {1}'.format(job_dict['recepmod_fn'], frame_dir)
      subprocess.call(cmd_st, shell = True, cwd = prefix, stdout=subprocess.PIPE,
                      stdin=subprocess.PIPE)

    job_dict['docked_fn'] = '{0}_level{1}_docked_{2}.sdf'.format(job_dict['docker'], job_dict['docker_level'], job_dict['compoundstate_id'])
    # Copy the docked file to the frame dir
    cmd_st = 'cp {0} {1}'.format(job_dict['docked_fn'], frame_dir)
    subprocess.call(cmd_st, shell = True, cwd = prefix, stdout=subprocess.PIPE,
                      stdin=subprocess.PIPE)
    job_dict['dest_dir']        = os.path.join(local_paths['resource_prefix'],
                                        job_dict['dbname'])
    job_dir = os.path.join(job_dict['dest_dir'], 'job'+str(job_dict['job_id']))
    if not os.path.exists(job_dir):
      os.makedirs(job_dir)

    
    #save everything to job_dir
    cmd_st = 'cp * {0}'.format(job_dir)
    subprocess.call(cmd_st, shell = True, cwd = prefix, stdout=subprocess.PIPE,
                    stdin=subprocess.PIPE)
    #with open(os.path.join(job_dir, 'job_dict.txt'), 'w') as ofp:
    #  ofp.write(str(d))
    score_fn = os.path.join(job_dir, 'score.txt')

    results = lddb_parsers.parse_fred2(score_fn)
    
    st_out = cStringIO.StringIO()
    for r in results:
      r = [job_dict['job_id']] + r
      r = map(str, r)
      st_out.write('\t'.join(r) + '\n')

    st_out.seek(0)
    cur = conn.cursor()
    cur.copy_from(st_out, 'Level1_FredScores', columns = ('job_id', 'cg4', 'steric', 'clash', 'pro_desolv', 'lig_desolv', 'lig_desolv_hb', 'cg4_hb'))
    #cur.copy_from(st_out, 'Level1_FredScores', columns = ('job_id', 'gs', 'gvdw', 'ges', 'inen'))
    
    cur.close()
    conn.commit()
    return bool(results), job_dir

#
#=====================================================================================



def test_fred(input_dict, local_paths):
  #print '\tTesting fred'
  FredDocker.preprocess(input_dict, local_paths)
  #print '\tPreprocessed'
  gen = FredDocker
  #print '\tMade gen'
  FredDocker.runFred(input_dict)
  print '\tran'
  #FredDocker.load(conn, prefix, input_dict, local_paths)
  #how to do conn?  

if __name__ == '__main__':
  '''
  input_dict = { "m_id": 1, 
      "gromacs_pdb2gmx": "/damsl/projects/molecules/software/MD/Gromacs/gromacs/gromacs-4.5.5/gromacs/bin/pdb2gmx", 
      "run_at": "molecules/data/Odorant_GPCR/", 
      "gromacs_grompp": "/damsl/projects/molecules/software/MD/Gromacs/gromacs/gromacs-4.5.5/gromacs/bin/grompp", 
      "local_prefix": "/damsl/projects/", 
      "job_id": 23, 
      "generator": "fred", 
      "charmm_bin": "/damsl/projects/molecules/software/MD/Charmm/c36b1/exec/gnu/charmm", 
      "resource_prefix": "/damsl/projects/", 
      "c_id": 12, 
      "r_id": 1393, 
      "recepmod_fred_dir": "molecules/data/Odorant_GPCR/receptor_models/olfr1393/oe_fred_high_gpcrm_model1", 
      "recepmod_fred_fn": "gpcrm_model1_olfr1393_grid.oeb.gz", 
      "compound_dir": "molecules/data/Odorant_GPCR/compounds", 
      "smi_fn": "ODL00000012_10022-28-3.smi", 
      "recepmod_dock6_fn": "gpcrm_model_r01_mod2_aligned.grid.bmp", 
      "omega2_bin": "/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/omega2", 
      "io_dir": "local_io", 
      "recepmod_dock6_dir": "molecules/data/Odorant_GPCR/receptor_models/olfr1393/dock6_finegrid_gpcrm_model1", 
      "gromacs_mdrun": "/damsl/projects/molecules/software/MD/Gromacs/gromacs/gromacs-4.5.5/gromacs/bin/mdrun", 
      "taut": False, 
      "namd_bin": "/damsl/projects/molecules/software/MD/NAMD/NAMD/NAMD_2.9_Linux-x86_64-multicore/namd2",
  #
  '''
  input_dict = {
       #'protein_pdb': '/damsl/projects/molecules/data/odorant_gpcr/membrane_bilayer_md/charmm-gui/gpcrm_model_r03_proa.pdb',
       'resource_prefix': '/damsl/projects',
       'io_dir':          'resource_io',
       'deployment_dir': 'd_00_1',
       'session_dir':     'dummy_session',
       'run_dir':         'run_2',
       'template_prefix': '/home/kwong23/assembla/damsl.mddb/templates',  
       #'template_prefix': '/home/nutanong/mddb/templates/',  
       'trj_dir':         '/damsl/projects/molecules/data/Odorant_GPCR/receptor_models/olfr1393/oe_fred_high_gpcrm_model6/bilayer_frames',
       'pdb_fn':          'protein.pdb',
       'smiles':          'C(C)(C)(C)C1CCC(=O)CC1',
       'trj_id':          '21',
       't':               '1',
       'model_id':        '7',
       'receptor_id':     '1393',
       'docker':          'fred',
       'docker_level':    '1',
       'compoundstate_id': '1',
       'dest_dir':        '/mddb/project/data/lddb_odorant2',
       'chimera':         '/damsl/projects/molecules/software/DD/Chimera/exec/bin/chimera',
       'sphgen_bin':           '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/sphgen',
       'molcharge_bin':           '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/molcharge',
       'dock6_bin':           '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/dock6',
       'dbname':          'lddb_odorant2',
       'sphgen':          '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/sphgen',
       'showbox':         '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/showbox',
       'grid':            '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/grid',
       'molcharge':       '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/molcharge',
       'omega2':          '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/omega2',
       'dock6':           '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/dock6',
       'sphere_selector': '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/sphere_selector',
       'omega2_bin':           '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/omega2',

      'fred_bin': '/damsl/projects/molecules/software/DD/OpenEye/centos/openeye/bin/fred', 
      'apopdb2receptor': '/damsl/projects/molecules/software/DD/OpenEye/openeye/arch/redhat-RHEL6-x64/oedocking/3.0.1/bin/apopdb2receptor'
      }


  try:

    print 'Begin test.'
    import resources
    #print 'Resources imported.'
    test_fred(input_dict, resources.LocalResource.get_paths())
    print 'Done test.'
  except:
    pass





