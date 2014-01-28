
from abc import ABCMeta, abstractmethod, abstractproperty
import os
import mddb_utils
#import generate_peptides as genpep
import datetime, time
import random
import math
import cStringIO
import generator
import subprocess
import md_pipelines
import shutil

try:
  import parser as mddb_parser
except:
  pass

class AbstractSimulator(generator.AbstractGenerator):
  __metaclass__ = ABCMeta

  @abstractmethod
  def run(self, input_params):
    raise NotImplementedError( "Should have implemented this" )

  @staticmethod
  def get_output_fns(job_dict):
    raise NotImplementedError( "Should have implemented this" )


  @staticmethod
  def get_sync_info(input_dict):
    raise NotImplementedError( "Should have implemented this" )

  def run_min(self, script_params):
    done_min = False
    count = 1
    start_from_anywhere = False
  
  
    start_from_anywhere = False
    phi_psi_array = script_params['phi_psi_array']
    if len(phi_psi_array) == 0:
      phi_psi_array = map(lambda x: (-math.pi + random.random()*math.pi*2),
                                  [None] * (len(script_params['protein_seq'])-1)*2)
      script_params['phi_psi_array'] = phi_psi_array
      start_from_anywhere = True
      
    

    """Run Minimization using Charmm by default"""
    new_job_dict = None
    while True:
      new_job_dict = genpep.generate_scripts(script_params)

      pdb_fn = new_job_dict['pdb_file_name']
      if os.path.isfile(pdb_fn):
        mddb_utils.run_cmd_with_retries("rm {0}".format(pdb_fn), 1000)

      minimization_in_fn = new_job_dict['minimization_in_fn']
      while True:
        try:
          charmm_bin   = script_params['charmm_bin']
          st =mddb_utils.run_cmd_with_file(charmm_bin,
                                 minimization_in_fn, minimization_in_fn+".log")
          break
        except Exception, e:
          print str(e)
          time.sleep(0.2)

      if not os.path.isfile(pdb_fn):
        time.sleep(0.5)
      else:
        done_min = os.path.getsize(pdb_fn) > 0

      if done_min:
        break

      neg_or_pos = [-1,1]

      print "Trial {0}: {1} {2}".format(count, script_params['phi_psi_array'], done_min)
      count = count + 1
      if count < 50 and not start_from_anywhere:
        delta_phipsi =  map(lambda x: neg_or_pos[random.randint(0,1)] *  (random.uniform(0,10)*math.pi/180),
                                    [None] * (len(script_params['protein_seq'])-1)*2)
        script_params['phi_psi_array'] = map(lambda a,b: a-b, script_params['phi_psi_array'], delta_phipsi)
      else:
        sign = neg_or_pos[random.randint(0,1)]
        script_params['phi_psi_array'] = map(lambda x: sign* (-math.pi + random.randint(4,15)*math.pi/24),
                                    [None] * (len(script_params['protein_seq'])-1)*2)

    print new_job_dict['simulation_in_fn']
    return dict(script_params.items() + new_job_dict.items())

  @staticmethod
  def reformat_velocity_dcd(vel_fn):
    with open(vel_fn, 'rb') as ifp:
      content = ifp.read()

    content = content.replace(b'VELD', b'CORD', 1)

    with open(vel_fn, 'wb') as ofp:
      ofp.write(content)

  def load(self, conn, result_dir, d, local_paths):
    d['prefix'] = result_dir
    out_fns = self.__class__.get_output_fns(d)

    trj_id = d['trj_id']
    e_id   = d['e_id'] 
    sim_fn = out_fns['sim_fn']
    trj_fn = out_fns['trj_fn']
    nrg_fn = out_fns['nrg_fn']
    vel_fn = out_fns.get('vel_fn')

    print "simulation file: ", sim_fn
    print "trajectory file: ", trj_fn
  
    ap_t_name = 'trj_{0}_ap'.format(trj_id)
    pf_t_name = 'trj_{0}_pf'.format(trj_id)
  
    step = 1
    l = mddb_parser.get_trjectory_length(sim_fn, trj_fn)
    chunk_size = 100000
    num_chunks = int(math.ceil(float(l)/chunk_size))
    cur = conn.cursor()
    for i in xrange(0, num_chunks):
      start = i * chunk_size
      stop  = min(start + chunk_size, l)
  
      print "trj_id: ", trj_id, "      range: ", start, "->", stop \
            ,"    total: ", l, "         time: ", datetime.datetime.now()
  
      crd_st_io = mddb_parser.parse_trajectory_range_to_buffer(trj_id, sim_fn, trj_fn, start, stop, step)

      if vel_fn and os.path.isfile(vel_fn):
        print "velocity file: ", vel_fn
        if os.path.splitext(vel_fn)[1] == '.dcd':
          self.__class__.reformat_velocity_dcd(vel_fn)

        vel_st_io = mddb_parser.parse_trajectory_range_to_buffer(trj_id, sim_fn, vel_fn, start, stop, step)
        if vel_st_io != None:
          cur.copy_from(vel_st_io, 'AtomVelocities', columns=('trj_id', 't', 'atom_id', 'x', 'y', 'z'))
  
      if crd_st_io == None:
        print "Empty string returned"
        return False
  
      st = "select create_NewAtomPositions_table(E'{0}', {1});".format(ap_t_name, trj_id)
      cur.execute(st)
      cur.copy_from(crd_st_io, ap_t_name, columns=('trj_id', 't', 'atom_id', 'x', 'y', 'z'))
  
      cur.execute('create index {0}Index on {0} (trj_id, t);'.format(ap_t_name))
      cur.execute("select compute_features ({0}, E'{1}', E'{2}')".format(e_id, ap_t_name, pf_t_name))
      cur.execute('insert into AtomPositions select * from {0}'.format(ap_t_name))
      cur.execute('insert into ConformationSpace select * from {0}'.format(pf_t_name))
  
      cur.execute("drop table if exists {0}".format(ap_t_name))
      cur.execute("drop table if exists {0}".format(pf_t_name))
  
  
    with open(nrg_fn, 'r') as ifp:
      table_width = len(ifp.readline().split())
      ifp.seek(0)
      if table_width == 3:
        cur.copy_from(ifp, 'Energies', columns=('trj_id', 't', 'total'))
      elif table_width == 16:
        cur.copy_from(ifp, 'Energies', columns=('trj_id', 't', 'bond', 'angle', 'dihed',
                                                'imprp', 'elect', 'vdw', 'boundary',
                                                'misc', 'kinetic', 'total', 'temp', 'potential',
                                                'total3', 'tempavg'))
    cur.close()
    conn.commit()
    return True
  


#=====================================================================================


class CharmmSimulator(AbstractSimulator):
  def run(self, input_params):
    trj_id = input_params['trj_id']
    new_dict = self.run_min(input_params)
    print input_params
    in_fn  = new_dict['simulation_in_fn']
    bin    = new_dict['charmm_bin']
    log_fn = in_fn + '.log'
    nrg_fn = log_fn + '.nrg'
    mddb_utils.run_cmd_with_file(bin, in_fn, log_fn)

    st_io = self.parse_energies(log_fn, trj_id)
    with open(nrg_fn, 'w') as ofp:
      ofp.write(st_io.getvalue())


  def parse_energies(self, log_fn, trj_id):
    print 'parse charmm energy value'

    #log_fn = d['log_fn']
    #trj_id = d['trj_id']

    output_st   = cStringIO.StringIO()
    t=0
    with open(log_fn, "r") as inp:
      for line in inp:
        if t == 0:
          t = 1
          continue
        if "DYNA>" in line[0:5]:
          output_line = "\t".join([str(trj_id), str(t-1)] + line.split()[3:4]) + "\n"
          output_st.write(output_line)
          t=t+1
    output_st.seek(0)
    return output_st

  @staticmethod
  def get_output_fns(job_dict):
    script_params = dict(job_dict)
    script_params['fns_only'] = True
    fn_dict = genpep.generate_scripts(script_params)
    out_dict = { 'sim_fn': fn_dict['psf_file_name']
                ,'pdb_fn': fn_dict['pdb_file_name']
                ,'log_fn': fn_dict['simulation_in_fn'] + '.log'
                ,'trj_fn': fn_dict['dcd_file_name']
                ,'vel_fn': fn_dict['vel_file_name']
                ,'nrg_fn': fn_dict['simulation_in_fn'] + '.log.nrg'
               }
    return out_dict

  @staticmethod
  def get_sync_info(input_dict):
    return []
#=====================================================================================


class NAMDSimulator(AbstractSimulator):
  def run(self, input_params):
    print input_params
    in_fn = input_params['simulation_in_fn']
    bin   = input_params['namd_bin']
    log_fn = input_params['simulation_in_fn'] + '.log'
    mddb_utils.run_cmd_with_env(bin + " " + in_fn, log_fn, None)
    ret_dict = \
      { 'trj_fn': input_params['out_file_name'] + ".dcd"
       ,'sim_fn': input_params['psf_file_name']
       ,'log_fn': log_fn
       ,'trj_id': input_params['trj_id']
      }

    st_io = generator.parse_energies(out_dict)
    with open(ret_dict['nrg_fn'], 'w') as ofp:
      ofp.write(st_io.getvalue())

    return ret_dict

  def parse_energies(self, d):
    log_fn = d['log_fn']
    trj_id = d['trj_id']
    output_st   = cStringIO.StringIO()
    t=0
    with open(log_fn, "r") as inp:
      for line in inp:
        if t == 0:
          t = 1
          continue
        if "ENERGY:" in line[0:7]:
          #print t-1, line
          output_line = "\t".join([str(trj_id), str(t)] + line.split()[2:]) + "\n"
          #print output_line
          output_st.write(output_line)
          t=t+1
    output_st.seek(0)
    return output_st

#=====================================================================================

class GromacsSimulator(AbstractSimulator):
  def gen_conversion_script(self,in_pdb_fn, out_pdb_fn, script_fn):

    st = """library(bio3d)\n""" +\
         """pdb <- read.pdb('{0}')\n""".format(in_pdb_fn) +\
         """new <- convert.pdb(pdb, type="amber")\n""" +\
         """write.pdb(new, file='{0}')\n""".format(out_pdb_fn)

    open(script_fn, 'w').write(st)

  def reformat_pdb1(self, pdb_fn):
    mddb_utils.run_cmd("cp {0} {0}.bak".format(pdb_fn))
    script_fn =  pdb_fn + '.R'
    self.gen_conversion_script(pdb_fn, pdb_fn + '2', script_fn)
    mddb_utils.run_cmd("Rscript " + script_fn)
    mddb_utils.run_cmd("mv {0}2 {0}".format(pdb_fn))


  def run(self, input_params):
    print input_params
    in_fn      = input_params['simulation_in_fn']
    pdb2gmx    = input_params['gromacs_pdb2gmx']
    grompp     = input_params['gromacs_grompp']
    mdrun      = input_params['gromacs_mdrun']

    gro_fn     = input_params['gro_file_name']
    tpr_fn     = input_params['tpr_file_name']
    xtc_fn     = input_params['xtc_file_name']
    trr_fn     = input_params['trr_file_name']
    pdb_fn     = input_params['pdb_file_name']
    mdp_fn     = input_params['simulation_in_fn']
    top_fn     = input_params['top1_file_name']

    self.reformat_pdb1(pdb_fn)
    mddb_utils.run_cmd(pdb2gmx + " -f " + pdb_fn + " -o " + gro_fn +\
                       " -p " + top_fn + "-ignh -ff amber99sb -water none")

    if not os.path.isfile(top_fn):
      top_fn = input_params['top2_file_name']

    mddb_utils.run_cmd(grompp + " -v -f " + mdp_fn + " -c " + gro_fn +\
                       " -p " + top_fn + " -o " + tpr_fn)

    tmp_xtc_fn = "/tmp" + xtc_fn;
    tmp_trr_fn = "/tmp" + trr_fn;

    try:
      os.makedirs(os.path.dirname(tmp_xtc_fn))
    except:
      pass

    log_fn = input_params['simulation_in_fn'] + '.log'
    sim_cmd = mdrun + " -s " + tpr_fn + " -x " + tmp_xtc_fn + " -o " + tmp_trr_fn + " > " + log_fn

    mddb_utils.run_cmd(sim_cmd)

    mddb_utils.run_cmd("cp {0} {1}".format(tmp_xtc_fn, xtc_fn))
    mddb_utils.run_cmd("cp {0} {1}".format(tmp_trr_fn, trr_fn))

    ret_dict = \
      { 'trj_fn': trr_fn
       ,'sim_fn': gro_fn
       ,'log_fn': log_fn
      }
    return ret_dict


#=====================================================================================
# Amber explicit solvent simulator 
class AmberSimulator(AbstractSimulator):
  
  # Run Amber simulation pipeline
  def run(self, input_params):
    pl = md_pipelines.AmberPipeline(input_params)
    pl.run()
    with open(os.path.join(input_params['output_prefix'], 'status.txt'), 'w') as ofp:
      ofp.write('done')

  # Copy data to a local running location and performs a light computation e.g., sed command
  @staticmethod
  def preprocess(d, local_paths):
    md_pipelines.AmberPipeline.preprocess(d, local_paths)

  # expected output files.. the deployment will wait for this one
  @staticmethod
  def get_output_fns(d):
    out_dict = { 'sta_fn': 'status.txt' }
    return out_dict

  # data locations to be synchronized to reduce data movement per job
  @staticmethod
  def get_sync_info(d):
    return []

  # load the data by copying the run/result directory to an appropriate location
  def load(self, conn, result_dir, d, local_paths):
    t = d.get('t')
    if t != None:
      dest_dir = os.path.join(d['dest_dir'], "trj{trj_id}_t{t}_ns{nstep_simulation}".format(**d))
    else:
      dest_dir = os.path.join(d['dest_dir'], "trj{0}".format(d['trj_id']))


    while True:
      try:
        if os.path.isdir(dest_dir):
          shutil.rmtree(dest_dir)

        shutil.copytree(result_dir, dest_dir)
        break
      except:
        pass

    return True,dest_dir

#
# Membrane simulator for which uses Charmm and NAMD minimizations and NAMD for production
#
class MembraneSimulator(AbstractSimulator):

  def run(self, input_params):
    input_params['protein_pdb'] = 'protein.pdb'
    m = md_pipelines.MembranePipeline(input_params)
    m.run()
    output_prefix   = os.path.join(input_params['resource_prefix'],
                                   input_params['io_dir'],
                                   input_params['session_dir'],
                                   input_params['deployment_dir'],
                                   input_params['run_dir'])


    with open(os.path.join(output_prefix, 'status.txt'), 'w') as ofp:
      ofp.write('done')

  @staticmethod
  def preprocess(d, local_paths):
    md_pipelines.MembranePipeline.preprocess(d, local_paths)

  # expected output files.. the deployment will wait for this one
  @staticmethod
  def get_output_fns(d):
    out_dict = { 'sta_fn': 'status.txt' }
    return out_dict

  # data locations to be synchronized to reduce data movement per job
  @staticmethod
  def get_sync_info(d):
    return []

  # load the data by copying the run/result directory to an appropriate location
  def load(self, conn, result_dir, d, local_paths):
    recp_dir = os.path.join(d['dest_dir'], "membrane_runs/r{receptor_id}/".format(**d))

    if not os.path.exists(recp_dir):
      os.makedirs(recp_dir)

    model_dir = os.path.join(recp_dir,str(d['model_id']))
    print 'model_dir:', model_dir
    while True:
      try:
        if os.path.isdir(model_dir):
          shutil.rmtree(model_dir)

        shutil.copytree(result_dir, model_dir)
        break
      except:
        pass

    return True,model_dir

#
# 
#
#
#  @staticmethod
#  def get_output_fns(d):
#    return {  'sim_fn': 'step5_assembly.xplor_ext.psf'
#             ,'pdb_fn': None
#             ,'log_fn': 'namd/step7.1_production.log'
#             ,'trj_fn': 'namd/step7.1_production.dcd'
#             ,'vel_fn': None
#             ,'nrg_fn': None
#             ,'sta_fn': 'status.txt'
#           }
#
#  @staticmethod
#  def get_sync_info(d):
#    return []
#
#  def load(self, conn, result_dir, d, local_paths):
#    dest_dir = os.path.join(d['dest_dir'], "recp{receptor_id}_model{model_id}".format(**d))
#
#    if os.path.isdir(dest_dir):
#      shutil.rmtree(dest_dir)
#
#    shutil.copytree(result_dir, dest_dir)
#    return True,dest_dir
#
#


def TestAmberSim():
    d = \
      { 'output_directory':'/damsl/projects/resource_io/dummy_session/trj022_t1/',
        'simulation_in_fn':'/damsl/projects/resource_io/dummy_session/trj022_t1/amber_dyn_2e',
        'amber_bin': '/damsl/projects/molecules/software/MD/Amber/amber12/bin/pmemd.cuda_SPFP',
        'prm_file_name': 'protein.prmtop',
        'start_cond_fn': 'protein_pmdmd_cuda_tip4p_Ew_dyn1.rst',
        'end_cond_fn': 'protein_pmdmd_cuda_tip4p_Ew_dyn2.rst',
        'trj_fn': 'protein_pmdmd_cuda_tip4p_Ew_dyn2.ncdf'
      }
    mysim = AmberSimulator()

    return mysim.run(d)



if __name__ == '__main__':

  l = [CharmmSimulator]

  for c in l:
    print str(c),':'
    print '    Subclass:', issubclass(c, AbstractSimulator)
    print '    Instance:', isinstance(c(), AbstractSimulator)
   

  TestAmberSim() 
