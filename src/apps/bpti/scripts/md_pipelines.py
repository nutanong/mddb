from abc import ABCMeta, abstractmethod, abstractproperty
import os
import cStringIO
import pipeline
import mddb_utils
import subprocess
import time
import sys
import chunkIt
import math
import numpy

try:
  import re
  import MDAnalysis
except:
  pass


class MembranePipeline(pipeline.Pipeline):
  __metaclass__ = ABCMeta

  template_dir = 'membrane_min'
  def __init__(self, d):
    self.param_dict = d
    self.stages = self.param_dict.get('stages')
    d['template_dir'] = self.template_dir 

    # list of input files for dict sub
    self.inp_fn_list = [ 'toppar.str', 'step2_area.str', 'step2.1_pore.str'
                        , 'step3_nlipids_lower.prm', 'step2.1_pore_watbox.str'
                        , 'step3_nlipids_upper.prm', 'step3_packing_pol.str'
                        , 'crystal_image.str', 'step2.1_porewat.crd'
                        , 'step4.4_pore.str',  'step4_components.str'
                        , 'step5_assembly.str', 'checkfft.py'
                        , 'membrane_lipid_restraint.str'
                        , 'membrane_lipid_restraint2.str'
                        , 'step1_pdbreader.inp'
                        , 'step2_orient.inp'
                        , 'step2.1_pore.inp'
                        , 'step3_size.inp', 'step3_packing.inp'
                        , 'step4_lipid.inp', 'step4.2_waterbox.inp', 'step4.3_ion.inp', 'step4.4_pore.inp'
                        , 'step5_assembly.inp'
                        , 'membrane_lipid_restraint.namd.str'
                        , 'membrane_lipid_restraint2.namd.str'
                        , 'step6.1_equilibration.inp'
                        , 'step6.2_equilibration.inp'
                        , 'step6.3_equilibration.inp'
                        , 'step6.4_equilibration.inp'
                        , 'step6.5_equilibration.inp'
                        , 'step6.6_equilibration.inp'
                       ] 

    # aux files to be copied without sub
    self.aux_fn_list = [ 'namd/membrane_lipid_restraint.namd.col' ]

    # if output is not specified the default output extension is appended to the input file name
    self.default_out_ext    = '.pdb'

    # stages defined as a two-level nested lists
    # first substage of each stage acts as a checkpoint
    # for example, if 'charmm < step3_packing.inp' fails the execution will restart at 'charmm < step3_size.inp'
    self.stage_list  = [[{'cmd': 'charmm',   'si': 'step1_pdbreader.inp', 'a': 'toppar.str'}],
                        [{'cmd': 'charmm',   'si': 'step2_orient.inp'}, 
                         {'cmd': 'charmm',   'si': 'step2.1_pore.inp', 'a': 'step2.1_pore_watbox.str'} ],
                        [{'cmd': 'charmm',   'si': 'step3_size.inp', 'o': ['step3_size.str']}, 
                         {'cmd': 'charmm',   'si': 'step3_packing.inp'} ],
                        [{'cmd': 'charmm',   'si': 'step4_lipid.inp', 'o': ['step4_lipid.psf']}, 
                         {'cmd': 'charmm',   'si': 'step4.2_waterbox.inp', 'o': ['step4.2_waterbox.crd']},
                         {'cmd': 'charmm',   'si': 'step4.3_ion.inp'}, 
                         {'cmd': 'charmm',   'si': 'step4.4_pore.inp'} ],
                        [{'cmd': 'charmm',   'si': 'step5_assembly.inp'}],
                        [{'cmd': 'namd_mpi', 'i':  'step6.1_equilibration.inp', 'd': 'namd', 'o': ['step6.1_equilibration.dcd']}, 
                         {'cmd': 'namd_mpi', 'i':  'step6.2_equilibration.inp', 'd': 'namd', 'o': ['step6.2_equilibration.dcd']}, 
                         {'cmd': 'namd_mpi', 'i':  'step6.3_equilibration.inp', 'd': 'namd', 'o': ['step6.3_equilibration.dcd']},
                         {'cmd': 'namd_mpi', 'i':  'step6.4_equilibration.inp', 'd': 'namd', 'o': ['step6.4_equilibration.dcd']}, 
                         {'cmd': 'namd_mpi', 'i':  'step6.5_equilibration.inp', 'd': 'namd', 'o': ['step6.5_equilibration.dcd']}, 
                         {'cmd': 'namd_mpi', 'i':  'step6.6_equilibration.inp', 'd': 'namd', 'o': ['step6.6_equilibration.dcd']}],
                        [{'cmd': 'namd_mpi', 'i':  'step7.1_production.inp'   , 'd': 'namd', 'o': ['step7.1_production.dcd']}],
                       ]

  # Copy the protein pdb and perform file reformatting through pdb2pqr, add_sedid_PTEN, and several sed commands
  @staticmethod
  def preprocess(d, local_paths):
    receptor_dir = '/damsl/projects/molecules/data/Odorant_GPCR/receptor_models/'
    protein_pdb = os.path.join(receptor_dir,
                               "{0}{1}".format('olfr', d['receptor_id']),
                               d['pdb_fn'])

    output_prefix = d['output_prefix']
    if not os.path.exists(output_prefix):
      os.makedirs(output_prefix)

    cmd0 = 'python /damsl/projects/molecules/software/tools/pdb2pqr/pdb2pqr-1.8/pdb2pqr.py ' +\
          '--ffout=charmm --ff=charmm {0} {1}'.format(protein_pdb, 'protein0.pdb')

    cmd1 = 'python /damsl/projects/molecules/data/PTEN/Full/Charmm/ACE/charmm-gui-model1/add_segid_PTEN.py ' +\
           '-f protein0.pdb -o protein1.pdb -p pdb -s ' +\
           '/damsl/projects/molecules/data/PTEN/Full/Charmm/ACE/charmm-gui-model1/resids.dat'

    cmd2 = "sed 's/ ATOM/ATOM/' protein1.pdb | " +\
           "sed 's/CD1 ILE/CD  ILE/' | " +\
           "sed 's/1CBDISU/ CB CYS/' | " +\
           "sed 's/1SGDISU/ SG CYS/' > protein.pdb "

    print cmd0
    subprocess.call(cmd0, shell=True, cwd = output_prefix)
    print cmd1
    subprocess.call(cmd1, shell=True, cwd = output_prefix)
    print cmd2
    subprocess.call(cmd2, shell=True, cwd = output_prefix)

    os.remove(os.path.join(output_prefix, 'protein0.pdb'))
    os.remove(os.path.join(output_prefix, 'protein1.pdb'))


  # run substage method for two input file types:
  #   'i' - the input file is given as a command line argument
  #   'si' - the input file is given as a standard input
  def run_substage(self, output_prefix, substage):
    in_fn = substage.get('i') or substage.get('si')

    out_fns = substage.get('o')

    # certain commands which requires us to specify the number of cores has to be
    # reformatted at this point
    cmd   = self.param_dict[substage.get('cmd')]
    num_cores = self.param_dict['num_cores_per_node'] * self.param_dict['num_nodes']
    cmd   = cmd.format(nc = num_cores)

    log_fn = os.path.splitext(in_fn)[0] + '.log'
    subdir = substage.get('d') or ""
    workdir = os.path.join(output_prefix,subdir)

    print '==============================================================='
    print 'in_fn:', in_fn
    print 'output_prefix:', output_prefix
    print 'subdir:', subdir
    print 'workdir:', workdir
    print 'log_fn:', log_fn
    print 'num_cores:', num_cores


    with open(os.path.join(workdir, log_fn), "w") as debug_log:
      if substage.get('si'):
        print 'command:\n', cmd, " < ", in_fn
        with open(os.path.join(workdir, in_fn), 'r') as ifp:
          subprocess.call(cmd, shell=True, stdin=ifp, stdout=debug_log,
                          cwd = workdir,
                          stderr=subprocess.STDOUT)
      else:
        print 'command:\n', cmd, " ", in_fn
        subprocess.call(cmd + " " + in_fn, shell=True, stdout=debug_log,
                        cwd = workdir,
                        stderr=subprocess.STDOUT)
    print '==============================================================='
    sys.stdout.flush()
    sys.stderr.flush()


    time.sleep(0.2)
    try:
      subprocess.call("tail -n 60 {0}".format(log_fn), shell=True, cwd = workdir,
                        stderr=subprocess.STDOUT)
    except:
      pass
 
    print workdir
    print out_fns
    print map(lambda out_fn: os.path.isfile(os.path.join(workdir, out_fn)), out_fns)
    # the substage is considered alright when all output files are produced
    # returning True in this case
    return all(map(lambda out_fn: os.path.isfile(os.path.join(workdir, out_fn)), out_fns))


class AmberPipeline(pipeline.Pipeline):
  __metaclass__ = ABCMeta

  def __init__(self, d):
    self.param_dict = d
    self.inp_fn_list = ['leap.in', 'amber_min_1b', 'amber_min_2', 'amber_dyn_1', 'amber_dyn_2e', 'bpti_chk.rst']

    #echo '  51.2600000  51.2600000  51.2600000  90.0000000  90.0000000  90.0000000' >> bpti_trj3264_ts500.inpcrd

    self.stage_list = [[{'cmd': 'tleap',
                         'args': '-f leap.in',
                         'i': 'leap.in',
                         'o': []}],
                       [{'cmd': 'sander', 
                         'args': '-O -i amber_min_1b ' +\
                                    '-o min1bout ' +\
                                    '-p protein.prmtop ' +\
                                    '-c protein.rst ' +\
                                    '-r protein_min1.rst ' +\
                                    '-ref protein.rst',
                         'i': 'amber_min_1b',
                         'o': ['min1bout']}],
                       [{'cmd': 'sander',
                         'args': '-O -i amber_min_2 ' +\
                                    '-o amber_min_2_tip4p_Ew.out ' +\
                                    '-p protein.prmtop ' +\
                                    '-c protein_min1.rst ' +\
                                    '-r protein_min2.rst',
                         'i': 'amber_min_2',
                         'o': ['amber_min_2_tip4p_Ew.out']}],
                       [{'cmd': 'cuda_spfp', 
                         'args': '-O -i amber_dyn_1 ' +\
                                    '-o amber_pmemd_cuda_dyn_1_tip4p_Ew.out '+\
                                    '-p protein.prmtop '+\
                                    '-c protein_min2.rst ' +\
                                    '-r protein_pmdmd_cuda_tip4p_Ew_dyn1.rst ' +\
                                    '-x protein_pmemd_cuda_tip4p_Ew_dyn_1.mdcrd ' +\
                                    '-ref protein_min2.rst',
                         'i': 'amber_dyn_1',
                         'o': ['amber_pmemd_cuda_dyn_1_tip4p_Ew.out']}],
                       [{'func': self.reassign_vel,
                         'args': ('protein_pmdmd_cuda_tip4p_Ew_dyn1.rst', 'vel.txt', 'protein_pmdmd_cuda_tip4p_Ew_dyn1_vel.rst')}],
                       [{'cmd': 'cuda_spfp',
                         'args': '-O -i amber_dyn_2e ' +\
                                    '-o amber_pmemd_cuda_dyn_2_tip4p_Ew.out ' +
                                    '-p protein.prmtop ' +\
                                    '-c protein_pmdmd_cuda_tip4p_Ew_dyn1_vel.rst ' +\
                                    '-r protein_pmdmd_cuda_tip4p_Ew_dyn2.rst ' +\
                                    '-x protein_pmemd_cuda_tip4p_Ew_dyn_2.nc ',
                         'i': 'amber_dyn_2e',
                         'o': ['amber_pmemd_cuda_dyn_2_tip4p_Ew.out']}]
                      ]        


  def reassign_vel(self, output_prefix, rst1_fn, vel_fn, rst2_fn):
    with open(os.path.join(output_prefix, rst1_fn), 'r') as ifp:
      rst1 = ifp.readlines()
      vec_vel1 = []
      n = len(rst1)
      for l in rst1[n/2+1:n-1]:
        vec_vel1   = vec_vel1 + list(chunkIt.fixed_size(l.split(), 3))
      speeds1    = map(lambda v: math.sqrt(float(v[0])**2 + float(v[1])**2 + float(v[2])** 2), vec_vel1)



    with open(os.path.join(output_prefix, vel_fn), 'r') as ifp:
      vec_vel2 = eval(ifp.read())

      speeds2 = map(lambda v: math.sqrt(v[0]**2 + v[1]**2 + v[2] ** 2), vec_vel2)

    factor = numpy.mean(speeds1)/numpy.mean(speeds2)

    vec_vel2 = map(lambda v: map(lambda x: x*factor, v), vec_vel2)

    with open(os.path.join(output_prefix, rst2_fn), 'w') as ofp:
      n = len(rst1)
      ofp.write(''.join(rst1[0:n/2+1]))
      for i in range(0,len(vec_vel2)):
        ofp.write('{0:12.7f}{1:12.7f}{2:12.7f}'.format(*(map(float, vec_vel2[i]))))
        if i % 2 or i == len(vec_vel2) -1:
          ofp.write('\n')
      ofp.write(rst1[-1])

  @staticmethod 
  def calc_vel(dataDir, trj_id, t):
    myu = MDAnalysis.Universe(dataDir+'bpti_from_mae.pdb',dataDir+'bpti-all/bpti-all-'+trj_id+'.dcd')
    myu.trajectory[t]
    s = myu.selectAtoms("not name pseu")
    cs2 = s.atoms.coordinates()

    t = t-1
    if t == -1:
      trj_id = "{0:03d}".format(int(trj_id) -1)
      myu = MDAnalysis.Universe(dataDir+'bpti_from_mae.pdb',dataDir+'bpti-all/bpti-all-'+trj_id+'.dcd')

    myu.trajectory[t]
    s = myu.selectAtoms("not name pseu")
    cs1 = s.atoms.coordinates()

    vec_vel = []
    vec_pos = []
    for i in range(len(cs2)):
        vec_pos.append(cs2[i])
        vec_vel.append(cs2[i]-cs1[i]) #converting to amber unit

    vec_vel = map(lambda x: x.tolist(), vec_vel)
    return vec_vel
  
#  @staticmethod
#  def output_vel(out_fn, myvel):
#    with open(out_fn, 'w') as ofp:
#      for i in range(0,len(myvel)-1,2):
#        ofp.write('{0:12.7f}'.format(myvel[i][0])+\
#                  '{0:12.7f}'.format(myvel[i][1])+\
#                  '{0:12.7f}'.format(myvel[i][2])+\
#                  '{0:12.7f}'.format(myvel[i+1][0])+\
#                  '{0:12.7f}'.format(myvel[i+1][1])+\
#                  '{0:12.7f}'.format(myvel[i+1][2])+'\n')
#  
#      ofp.write('{0:12.7f}'.format(myvel[-1][0])+\
#                '{0:12.7f}'.format(myvel[-1][1])+\
#                '{0:12.7f}'.format(myvel[-1][2])+'\n')
  

  @staticmethod
  def output_vel(out_fn, myvel):
    with open(out_fn, 'w') as ofp:
      #print myvel
      #print '=========================================='
      #print str(myvel)
      ofp.write(str(myvel))


  
  @staticmethod
  def preprocess(d, local_paths):
    trj_id      = d['trj_id']
    t           = d['t']

    output_prefix = d['output_prefix']
    if d['source'] == 'deshaw':
      dataDir = '/damsl/mddb/data/deshaw/MD0_DATA/Public/outgoing/science2010/DESRES-Trajectory-bpti-all/'
      if not os.path.exists(output_prefix):
        os.makedirs(output_prefix)

      vec_vel = AmberPipeline.calc_vel(dataDir, trj_id, t)
      #print "vec_vel", vec_vel
      AmberPipeline.output_vel(os.path.join(output_prefix, 'vel.txt'), vec_vel)

      sys.stdout.flush()
      sys.stderr.flush()
      myu = MDAnalysis.Universe(dataDir+'bpti_from_mae.pdb',dataDir+'bpti-all/bpti-all-'+trj_id+'.dcd')
      npu = myu.selectAtoms('not name pseu')
      if t:
          myu.trajectory[t]

      npu.write(os.path.join(output_prefix, 'temp'+trj_id+'.pdb'))

      pfile = open(os.path.join(output_prefix, 'temp'+trj_id+'.pdb'),'r')
      nfile = open(os.path.join(output_prefix,'protein.pdb'),'w')

      for line in pfile:
          tm = re.search(r' [0-9][A-Z][A-Z] ',line)
          if not tm:
              if 'CL' in line:
                  nfile.write(line.replace('CL ','Cl-',2))
              elif 'CYS' in line:
                  nfile.write(line.replace('CYS','CYX'))
              else:
                  nfile.write(line)
      pfile.close()
      os.remove(os.path.join(output_prefix, 'temp'+trj_id+'.pdb'))
      nfile.close()


  def run_substage(self, output_prefix, substage):
    out_fns = substage.get('o')
    in_fn   = substage.get('i')


    func   = substage.get('func')
    if func != None:
      args   = substage.get('args')
      func(output_prefix,*args)
      return True
    else:
      cmd   = self.param_dict.get(substage.get('cmd'))
      args  = substage.get('args').format(**self.param_dict)
  
      print '==============================================================='
      print 'in_fn:', in_fn
      print 'output_prefix:', output_prefix
      print 'command:\n', cmd, " ", args
      print 'out_fns: ', ', '.join(out_fns)
      print '==============================================================='
  
      sys.stdout.flush()
  
      log_fn = (in_fn or 'temp') + '.log'
      with open(os.path.join(output_prefix, log_fn), "w") as debug_log:
  
        subprocess.call(cmd + ' ' + args, shell=True, cwd = output_prefix,
                        env={'AMBERHOME':os.getenv('AMBERHOME'),
                             'LD_LIBRARY_PATH': os.getenv('LD_LIBRARY_PATH') or "/user/lib",
                             'CUDA_VISIBLE_DEVICES': str(self.param_dict.get('device_id'))
                            },
                        stderr=subprocess.STDOUT, stdout=debug_log)
      time.sleep(0.1)
  
      out_fns.append(log_fn)
      for ofn in out_fns:
        print ofn, '------------------------------------------------------'
        subprocess.call("tail -n 40 {0}".format(ofn), shell=True, cwd = output_prefix,
                         stderr=subprocess.STDOUT)
  
  
      sys.stdout.flush()
      #return all(map(lambda out_fn: os.path.isfile(os.path.join(output_prefix, out_fn)), out_fns))
      return True


def testCharmmMinimizer():
  d = {'exp_solv_toppar_fn': 'step0_toppar_str.charmm.inp',
       'protein_pdb': '/damsl/projects/molecules/data/odorant_gpcr/membrane_bilayer_md/charmm-gui/gpcrm_model_r03_proa.pdb',
       'resource_prefix': '/damsl/projects',
       'io_dir':          'resource_io',
       'session_dir':     'dummy_session',
       'run_dir':         'run_1',
       'template_dir':    '/home/nutanong/mddb/templates/membrane_min',
       'aux_dirs':      ['toppar','lipid_lib','namd'],
       'charmm':        'charmm',
       'namd':          '/damsl/projects/molecules/software/MD/NAMD/NAMD/namd2',
       'namd_mpi':    '/home1/00288/tg455591/NAMD_2.8/NAMD_2.8_Linux-x86_64-ibverbs-Lonestar/charmrun +p108 ++verbose ++mpiexec ++remote-shell /home1/00288/tg455591/NAMD_scripts/mpiexec ++runscript tacc_affinity /home1/00288/tg455591/NAMD_2.8/NAMD_2.8_Linux-x86_64-MVAPICH-Intel-Lonestar/namd2',
       'charmm':      '/home1/01654/twoolf/c36a2_xlarge_dims_lonestar',

      }

  minimizer = CharmmMinimizer(d)
  print minimizer.states2qsub()


def testAmberMinimizer():
  d = {'exp_solv_toppar_fn': 'step0_toppar_str.charmm.inp',
       'protein_pdb': '/damsl/projects/molecules/data/odorant_gpcr/membrane_bilayer_md/charmm-gui/gpcrm_model_r03_proa.pdb',
       'resource_prefix': '/damsl/projects',
       'io_dir':          'resource_io',
       'session_dir':     'dummy_session',
       'template_dir':    '/home/jjeliaz1/mddb/templates/amber_min',
       'run_dir':         'run_2',
       'protein_name':    'bpti',
       'aux_dirs':      [],
       'sander':          '/damsl/projects/molecules/software/MD/Amber/amber12/bin/sander',
       'cuda_spfp':       '/damsl/projects/molecules/software/MD/Amber/amber12/bin/pmemd.cuda_SPFP',
       'tleap':           '/damsl/projects/molecules/software/MD/Amber/amber12/bin/tleap',
       'echo':            'echo',
      }

  minimizer = AmberMinimizer(d)
  minimizer.run()


if __name__ == '__main__':
  testCharmmMinimizer()




