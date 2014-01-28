from abc import ABCMeta, abstractmethod, abstractproperty
import os
import cStringIO
import pipeline
import mddb_utils
import subprocess
import time
import sys
import shutil
import uuid

class Dock6Pipeline(pipeline.Pipeline):
  __metaclass__ = ABCMeta

  template_dir = 'dock6'
  def __init__(self, d):
    print 'Pipeline init'
    print 'job_dict:', d
    self.param_dict = d
    #self.timeout = 7200 #2hr was not enough
    #self.timeout = 86400 #24hr to test -- 300/1639 failed from 24 hour limit
    self.timeout = 172800 
    d['template_dir'] = self.template_dir
    d['box_size'] = 8

    self.inp_fn_list = ['INSPH', 'step2_showbox.in', 'step3_grid.in', 'step4_dock6.in']
    self.aux_fn_list = ['align_all.py', 'get_centroid.py', 'place_helium.py', 
                        'writedms.py',  'remove_dmpc.py',  '1gzm.pdb', 'oe_license.txt']

    self.stage_list = [[{'cmd': 'chimera', 
                         'args': '--nogui --silent align_all.py'.format(**d),
                         'o': 'protein_aligned.pdb'},
                        {'cmd': 'molcharge',
			 'args': '-method gasteiger protein_aligned.pdb protein_aligned_charged.mol2'.format(**d),
                         'o': 'protein_aligned_charged.mol2'},
                        {'cmd': 'chimera', 
                         'args': '--nogui --script "get_centroid.py"'.format(**d),
                         'o': 'protein_aligned_charged_noH.mol2',
                         'so': 'com.txt'},
                        {'cmd': 'chimera', 
                         'args': '--nogui --silent --script "place_helium.py protein_aligned.pdb"'.format(**d),
                         'o': 'protein_aligned.mol2'},
                        {'cmd': 'chimera',
                         'args': '--nogui --silent protein_aligned_charged_noH.mol2 writedms.py'.format(**d),
                         'o': 'rec2.dms'},
                        {'cmd': 'python', 
                         'args': 'remove_dmpc.py'.format(**d),
                         'o': 'rec.dms'},
                        {'cmd': 'sphgen',
                         'args': '-i INSPH -o OUTSPH',
                         'o': ['OUTSPH', 'temp1.ms', 'temp2.sph', 'temp3.atc', 'temp.log', 'rec.sph']},
                        {'cmd': 'sphere_selector',
                         'args': 'rec.sph protein_aligned.mol2 3',
                         'o': 'selected_spheres.sph' }],
                       [{'cmd': 'showbox',
                         'args': '< step2_showbox.in',
                         'o': ['step2_box.pdb']}],
                       [{'cmd': 'grid',
                         'args': '-i step3_grid.in -o grid.out',
                         'o': 'grid.out'}],
                       #start here at step4 if grid is done already
                       [{'cmd': 'omega2',
                         'args': '-in ligand.smi -out ligand.mol2',
                         'o': 'ligand.mol2'},
			{'cmd': 'dock6', 
                         'args': '-i step4_dock6.in -o step4_dock6.out'.format(**d),
                         'o': 'step4_dock6.out' }],
                     
                      ]
    if d.get('grid_fn'):
      self.stage_list = self.stage_list[3:]                     


  @staticmethod
  def preprocess(d, local_paths):
    d['dest_dir']        = os.path.join(local_paths['resource_prefix'],
                                        d['dbname'])
    frame_dir = os.path.join(d['dest_dir'], "recp{receptor_id}/model{model_id}/trj{trj_id}/frame{t}".format(**d))
    #frame_dir     = os.path.join(recpmodel_dir, "frame{t}".format(**d))
    d['frame_dir'] = frame_dir
    if not os.path.exists(frame_dir):
      os.makedirs(frame_dir)

    frame_fn = os.path.join(frame_dir, 'protein.pdb')
    #copy pdb to the frame_dir
    if not os.path.exists(frame_fn):
      print 'trj_dir:', d['trj_dir']
      print 'pdb_fn:', d['pdb_fn']
      shutil.copy2(os.path.join(d['trj_dir'], d['pdb_fn']), frame_fn)

    output_prefix = d['output_prefix']  
    #create the run_dir
    if not os.path.exists(output_prefix):
      os.makedirs(output_prefix)

    #copy all previous frame_dir files to the run_dir
    cmd = 'cp {0}/* {1}/'.format(frame_dir, output_prefix)
    print cmd
    subprocess.call(cmd, shell=True)
    #write the smile file
    smile_fn = 'ligand.smi'
    with open(os.path.join(output_prefix, smile_fn), "w") as smi:
      smi.write(d['smiles'])
 
    #check if we have done the grid calculation already, if so, we can skip to step4
    grid_fn = 'grid.out'
    selected_spheres_fn = 'selected_spheres.sph'
    boxpdb_fn = 'step2_box.pdb'
    bmp_fn = 'step3_grid_score.bmp'
    nrg_fn = 'step3_grid_score.nrg'
    if os.path.exists(os.path.join(output_prefix, grid_fn)) and \
       os.path.exists(os.path.join(output_prefix, selected_spheres_fn)) and \
       os.path.exists(os.path.join(output_prefix, bmp_fn)) and \
       os.path.exists(os.path.join(output_prefix, nrg_fn)) and \
       os.path.exists(os.path.join(output_prefix, boxpdb_fn)):
      d['grid_fn'] = grid_fn
      print 'grid files exist'

   

  def run_substage(self, output_prefix, substage):
    in_fn   = substage.get('i') or substage.get('si')
    out_fns = substage.get('o') if isinstance(substage.get('o'), list) else [substage.get('o')]
    cmd     = self.param_dict.get(substage.get('cmd')) or substage.get('cmd')

    log_fn  = substage.get('so') or ((in_fn or 'temp') + '.log')
    subdir  = substage.get('d') or ""
    workdir = os.path.join(output_prefix,subdir)

    args    = substage.get('args')

    print '==============================================================='
    print 'in_fn:', in_fn
    print 'output_prefix:', output_prefix
    print 'workdir:', workdir
    print 'log_fn:', log_fn
    sys.stdout.flush()
    sys.stderr.flush()

    try:
      subprocess.call('rm ' + ' '.join(out_fns), shell=True,
                  cwd = workdir,
                  stderr=subprocess.STDOUT)
    except:
      pass

    with open(os.path.join(workdir, log_fn), 'w') as ofp:
      print 'command:\n', cmd, " ", args
      subprocess.call(cmd + " " + args, shell=True, stdout=ofp,
                      cwd = workdir,
                      stderr=subprocess.STDOUT)
    print '==============================================================='
    sys.stdout.flush()
    sys.stderr.flush()

    time.sleep(0.1)
    #subprocess.call("tail -n 60 {0}".format(log_fn), shell=True, cwd = output_prefix,
    #                  stderr=subprocess.STDOUT)

    return True



def test():
  d = {'protein_pdb': '/damsl/projects/molecules/data/odorant_gpcr/membrane_bilayer_md/charmm-gui/gpcrm_model_r03_proa.pdb',
       'resource_prefix': '/damsl/projects',
       'io_dir':          'resource_io',
       'session_dir':     'dummy_session',
       'run_dir':         'run_2',
       'template_prefix': '/home/kwong23/assembla/damsl.mddb/templates',  
       #'template_prefix': '/home/nutanong/mddb/templates/',  
       'trj_dir':         '/damsl/projects/molecules/data/Odorant_GPCR/receptor_models/olfr1393/oe_fred_high_gpcrm_model6/bilayer_frames',
       'pdb_fn':          'rosetta_modeller_model6_frame_0.pdb',
       'smiles':          'C(C)(C)(C)C1CCC(=O)CC1',
       'trj_id':          '21',
       't':               '1',
       'model_id':        '7',
       'receptor_id':     '1393',
       'dest_dir':        '/mddb/project/data/lddb_odorant2',
       'chimera':         '/damsl/projects/molecules/software/DD/Chimera/exec/bin/chimera',
       'sphgen_bin':           '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/sphgen',
       'molcharge_bin':           '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/molcharge',
       'omega2_bin':           '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/omega2',
       'dock6_bin':           '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/dock6',
       'dbname':          'lddb_odorant2',
       'sphgen':          '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/sphgen',
       'showbox':         '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/showbox',
       'grid':            '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/grid',
       'molcharge':       '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/molcharge',
       'omega2':          '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/omega2',
       'dock6':           '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/dock6',
       'sphere_selector': '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/sphere_selector',


      }

  try:
    import resources

    Dock6Pipeline.preprocess(d, resources.LocalResource.get_paths())

    pipeline = Dock6Pipeline(d)
    pipeline.run()
  except:
    pass


if __name__ == '__main__':
  test()



