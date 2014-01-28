# The resources module contains resource classes where each class
# provides an appropriate environment in which a data generator
# runs. For examples, 
#   - Environment variables PATH, PYTHONPATH, etc.
#   - Executables for Charmm, NAMD, Amber
#
# This removes machine dependency from data generation and makes
# the same data generation code runnable on different resources.

from abc import ABCMeta, abstractmethod, abstractproperty

import dockers
import simulators
import os
import platform
from optparse import OptionParser
import json
import subprocess
import uuid
import time
import threading
import random
import datetime
import mddb_utils
from multiprocessing import Process
import impsolv_simulators
import chunkIt
import cStringIO
import shutil
import sys


# A base resource class provides implementations for resource handling methods such as
# job deployment, data loading, data cleanup.

class Resource(object):

  __metaclass__ = ABCMeta

  res_configs_fn    = 'res_configs.txt'
  io_dir            = 'resource_io'
  local_prefix      = '/damsl/projects/mddb'
  script_dirs       = ['mddb/scheduler/','mddb/templates/']
  gateway_host      = None
  generator_options = { 'charmm':   simulators.CharmmSimulator,
                        'namd':     simulators.NAMDSimulator,
                        'gromacs':  simulators.GromacsSimulator,
                        'fred':     dockers.FredDocker,
                        'amber':    simulators.AmberSimulator,
                        'membrane': simulators.MembraneSimulator,
                        'dock6':    dockers.Dock6Docker,
                        'impsolv':  impsolv_simulators.ImpSolvSimulator,
                      }

  # Initialize resource configurations, e.g., #nodes per deployments, the name of the
  # job queue (if PBS), etc.
  def __init__(self, user, res_config_name, **kwargs):

    self.deployment_id = 1

    self.res_config_name     = res_config_name
    self.res_configs         = self.get_res_config(res_config_name, **kwargs)
    map(lambda x: setattr(self.__class__, *x), self.res_configs.items())

    try:
      gateway_host = self.gateway_host or 'localhost'
    except:
      gateway_host = 'localhost'

    self.user                = user
    self.process_pool        = []

    # a file name for keeping track of the last output sync (rsync) time
    # using it's modified time 
    self.sync_info_fn        = os.path.join('/tmp', gateway_host + '_' + str(uuid.uuid4()))
    # open and close the file to change the modified time
    open(self.sync_info_fn, 'w').close()

  # deployment id counter
  def get_next_deployment_id(self):
    deployment_id = self.deployment_id
    self.deployment_id = self.deployment_id + 1
    return deployment_id 

  # convert the sync period in secs into datetime.timedelta
  def get_sync_period(self):
    return datetime.timedelta(0, self.sync_period)


  # take an integer i; convert it to a string of a 6-digit hex; put a dash '-' every two characters
  # and put a prefix in front.
  # used to convert a deployment id and a run id into a directory name
  @staticmethod
  def int2dirname(dir_prefix, i):
    s = '{0:06x}'.format(i)
    return dir_prefix + '-'.join(s[i:i+2] for i in range(0, len(s), 2))

  @staticmethod
  def get_deployment_name(deployment_id):
    return Resource.int2dirname('d_', deployment_id)


  # get the resource config dictionary by modifying the default resource config dictionary
  # to get appropriate settings for MPI, GPU, normal resource configurations.
  def get_res_config(self, config_name, **kwargs):
    res_config = self.__class__.res_configs.get(config_name)
    if not res_config:
      return dict(self.default_config.items() + kwargs.items())
    return dict(self.default_config.items() + res_config.items() + kwargs.items())

  # open and close the file to change the modified time
  def set_last_sync(self):
    open(self.sync_info_fn, 'w').close()

  # read the modified time as a datetime object
  def get_last_sync(self):
    t = datetime.datetime.fromtimestamp( os.path.getmtime(self.sync_info_fn))
    return t

  # check for completed deployments to free up the deployment process pool
  # to allow more deployments. Returns true if the pool size is smaller than
  # the number of allowed deployments (e.g., the #queue entries limit for PBS)
  def check_deployments(self):
    finished_processes = [p for p in self.process_pool if not p.is_alive()]
    self.process_pool  = [p for p in self.process_pool if p not in finished_processes]

    for p in finished_processes:
      p.join()


    return len(self.process_pool) < self.num_deployments



  @staticmethod
  def get_environ():
    raise NotImplementedError( "Should have implemented this" )

  @staticmethod
  def print_environ(d):
    return " ".join(["{0}={1}".format(k,":".join(v)) for k,v in d.items()])

  # Starts a deployment process and insert it into the deployment process pool
  def deploy(self, session_dir, input_data, param_dict):
    deployment_id = self.get_next_deployment_id()
    self.prepare_job_dicts(session_dir, deployment_id, input_data, param_dict)
    p = Process(target = self.deploy_and_wait, args = (session_dir, deployment_id, input_data, param_dict))
    p.start()
    self.process_pool.append(p)

  # Add job running related attributes into each job dictionary
  # mostly the directory structure of where the job executions take place
  # on the resource. This method is resource invarient. There is no need
  # to overide/reimplement this definition.
  def prepare_job_dicts(self, session_dir, deployment_id, job_dicts, param_dict):
    st_io = cStringIO.StringIO()
    for i,d in enumerate(job_dicts):
      print d
      d.update(d.pop('data'))
      run_dir = self.int2dirname('r_', i)
      deployment_dir = self.get_deployment_name(deployment_id)

      d['run_dir'] = run_dir
      d['deployment_dir'] = deployment_dir
      d['deployment_id'] = deployment_id
      d['session_dir'] = session_dir
      d['user'] = self.user

      d['dbname']  = param_dict['dbname']

      local_paths = LocalResource.get_paths()

      d['output_prefix']   = os.path.join(local_paths['resource_prefix'],
                                          local_paths['io_dir'],
                                          session_dir,
                                          deployment_dir,
                                          run_dir)

      d['dest_dir']        = os.path.join(local_paths['resource_prefix'],
                                          d['dbname'])

      ts_str = "'{0}'".format(datetime.datetime.now())
      st_io.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(d['jq_entry_id'],
                                                self.worker_id,
                                                session_dir,
                                                deployment_id,
                                                i,
                                                ts_str
                                               ))

    st_io.seek(0)
    conn = mddb_utils.get_dbconn(param_dict['dbname'],
                                 param_dict['dbuser'],
                                 param_dict['dbhost'],
                                 param_dict['dbpass'])

    cur = conn.cursor()
    cur.copy_from(st_io, 'JobExecutionRecords', sep= '\t')
    conn.commit()
    conn.close()
    st_io.close()

  # Job deployment method (for remote resources)
  #   (1) iterate through a list of job dictionaries (job_dict_list) to preprocess each job
  #       (copying data to appropriate locations) and to find out expected
  #       output files.
  #   (2) Submit the job_dict_list to the associated remote resource
  #   (3) Wait until all expected output files are returend through rsync
  #   (4) Load the data to the database/file system
  #   (5) Clean up the remote and local job running directories

  def deploy_and_wait(self, session_dir, deployment_id, job_dict_list, param_dict):
    deployment_dir = self.get_deployment_name(deployment_id)
    output_data = []
    for d in job_dict_list:
      gen_name  = d.get('generator') or d.get('simulator') or d.get('docker')
      gen_class = Resource.generator_options[gen_name]
      gen_class.preprocess(d, LocalResource.get_paths())
      output_data.append((d['run_dir'],gen_class.get_output_fns(d)))
      d['user'] = self.user
    print 'job_dict_list:', len(job_dict_list), ' jobs'
    wait_list     = []
    for run_dir,d in output_data:
      for v in d.values():
        if v is not None:
          wait_list.append((run_dir, v))

    path_dict     = self.__class__.get_paths()
    remote_prefix = path_dict['resource_prefix']
    sync_dir      = os.path.join(Resource.io_dir, session_dir)
    try:
      os.makedirs(os.path.join(Resource.local_prefix,  sync_dir))
    except:
      pass

    self.submit(session_dir, deployment_id, job_dict_list)

    #print 'output files: ', wait_list
    while wait_list:
      time.sleep(60)

      ts = datetime.datetime.now()
      t_diff = ts - self.get_last_sync()
      if t_diff > self.get_sync_period():
        self.set_last_sync()
        print ts
        self.__class__.sync_output_dirs(self.__class__.gateway_host, [(sync_dir, None)], Resource.local_prefix, remote_prefix)
        self.set_last_sync()
        #print datetime.datetime.now()
        #print Resource.local_prefix
        #print sync_dir

      wait_list = [(run_dir,fn) for (run_dir,fn) in wait_list
                   if not os.path.isfile(os.path.join(Resource.local_prefix, sync_dir, deployment_dir, run_dir, fn))]
      if wait_list:
        print deployment_dir, 'is still waiting for', zip(*wait_list)[0]
      else:
        print deployment_dir, 'completed'

    conn = mddb_utils.get_dbconn(param_dict['dbname'],
                                 param_dict['dbuser'],
                                 param_dict['dbhost'],
                                 param_dict['dbpass'])


    self.load(conn, job_dict_list, LocalResource.get_paths())
    self.cleanup(job_dict_list, LocalResource.get_paths())

    conn.close()


  @staticmethod
  def get_paths():
    raise NotImplementedError( "Should have implemented this" )

  # Job loading method:
  #   (1) Iterate through a list of job dictionaries to call the appropriate load function
  #       according to the generator
  @staticmethod
  def load(conn, job_data, local_paths):

    if not isinstance(job_data, list):
      job_data = [job_data]

    st_io = cStringIO.StringIO()
    for d in job_data:
      print 'loading: ', d['jq_entry_id']
      result_dir = os.path.join(Resource.local_prefix, Resource.io_dir, 
                                d['session_dir'], d['deployment_dir'], d['run_dir'])
      print 'from: ', result_dir
      gen_name   = d.get('generator') or d.get('simulator') or d.get('docker')
      generator  = Resource.generator_options[gen_name]()

      try:
        ret = generator.load(conn, result_dir, d, local_paths)
      except:
        ret = (False, '')

      try:
        alright,path = ret
      except:
        alright = True
        path = ''


      ts_str = "'{0}'".format(datetime.datetime.now())
      st_io.write('{0}\t{1}\t{2}\t{3}\n'.format(d['jq_entry_id'],alright, path, ts_str))
      sys.stdout.flush()

    st_io.seek(0)
    cur = conn.cursor()
    cur.copy_from(st_io, 'JobCompletionRecords', sep= '\t')
    cur.close()
    conn.commit()
    st_io.close()

  # Data Clean Up method
  #   (1) remove the job running directory from the file system
  #   (2) if remote, ssh into the remote resource and remove the job running directory there as well

  def cleanup(self, job_data, local_paths):
    path_dict = self.__class__.get_paths()

    if not isinstance(job_data, list):
      job_data = [job_data]

    l = []
    for d in job_data:
      session_dir = d['session_dir']
      result_dir = os.path.join(Resource.local_prefix, Resource.io_dir, session_dir, d['deployment_dir'], d['run_dir'])
      try:
        shutil.rmtree(result_dir)
      except Exception as e:
        print type(e), e


      result_dir_remote = os.path.join(path_dict['resource_prefix'], 
                                       path_dict['io_dir'],
                                       session_dir,
                                       d['deployment_dir'],
                                       d['run_dir'])
      l.append(result_dir_remote)                                       
    if self.gateway_host:
      Resource.remove_remote_dirs(self.gateway_host, l)
      

  # Probing for an empty process slot ensuring that at any given point in time
  # the device id is unique which is critical for CUDA-based jobs
  @staticmethod
  def wait_avail_device(device_pool):
    while True:
      for device_id, proc in device_pool.items():
        if proc == None:
          return device_id
        elif not proc.is_alive():
          proc.join()
          return device_id
      time.sleep(2)
  

  # Execute n jobs in parallel where n = job_concurrency
  # If the number of jobs is greater than n, then some jobs
  # will have to wait for an empty slot. 
  def execute_jobs_parallel(self, job_data, path_dict):
    device_pool = dict([(i,None) for i in range(self.job_concurrency)])

    for d in job_data:
      device_id = Resource.wait_avail_device(device_pool)
  
      p = Process(target=self.execute_job, args=(d, path_dict, device_id))
      p.start()
      device_pool[device_id] = p

    for p in device_pool.values():
      try:
        p.join()
      except:
        pass

    return


  # execute jobs sequentially
  def execute_job(self, job_data, path_dict, device_id=0):
    
    if not isinstance(job_data, list):
      job_data = [job_data]

    out_data = []
    for input_dict in job_data:
      try:
        num_cores_per_node = self.num_cores_per_node
      except:
        num_cores_per_node = 1

      input_dict['num_cores_per_node'] = num_cores_per_node
      input_dict['device_id'] = device_id
      gen_name  = input_dict.get('generator') or input_dict.get('simulator') or input_dict.get('docker')
      generator = Resource.generator_options[gen_name]()

      run_dict = dict(input_dict.items() + path_dict.items() + self.res_configs.items())

      # specify the joblog file name
      output_prefix = os.path.join(run_dict['resource_prefix'], 
                                   run_dict['io_dir'], 
                                   run_dict['session_dir'], 
                                   run_dict['deployment_dir'],
                                   run_dict['run_dir'])

      joblog_fn = os.path.join(output_prefix, 'joblog.txt')
      jobdict_fn = os.path.join(output_prefix, 'jobdict.txt')
      with open(os.path.join(output_prefix, jobdict_fn), 'w') as ofp:
        ofp.write(str(input_dict))



      # save the original values for stdout and stderr 
      actual_stdout = sys.stdout
      actual_stderr = sys.stderr
      with open(joblog_fn,'w') as ofp:
        if self.__class__ != LocalResource:
          # redirect the stout and stderr to the joblog
          sys.stdout = ofp
          sys.stderr = ofp

        out_dict = generator.run(run_dict)

        if self.__class__ != LocalResource:
          # reset stdout and stderr to original values
          sys.stdout = actual_stdout
          sys.stderr = actual_stderr

      out_data.append(out_dict)


    if len(out_data) == 1:
      out_data = out_data[0]

    return out_data

  @staticmethod
  def sync_scripts(remote_host, remote_prefix, user, script_dirs):
    for script_dir in script_dirs:
      local_dir  = "{0}/{1}/{2}/".format(os.environ['HOME'], user, script_dir)
      remote_dir = "{0}/{1}/{2}/".format(remote_prefix, 'scripts', script_dir)
      Resource.sync_input(local_dir, remote_host, remote_dir, '*')

  @staticmethod
  def sync_input_dirs(remote_host, sync_list, local_prefix, remote_prefix):
    for dir, pattern in sync_list:
      local_dir  = "{0}/{1}/".format(local_prefix, dir)
      remote_dir = "{0}/{1}/".format(remote_prefix, dir)
      Resource.sync_input(local_dir, remote_host, remote_dir, pattern)

  @staticmethod
  def sync_output_dirs(remote_host, sync_list, local_prefix, remote_prefix):
    for dir, pattern in sync_list:
      local_dir  = "{0}/{1}/".format(local_prefix, dir)
      remote_dir = "{0}/{1}/".format(remote_prefix, dir)
      Resource.sync_output(local_dir, remote_host, remote_dir, pattern)

  @staticmethod
  def remove_remote_dirs(remote_host, remote_dirs):
    cmd_st = "ssh {0} rm -rf {1}".format(remote_host, ' '.join(remote_dirs))
    print cmd_st
    subprocess.Popen(cmd_st, shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()

  @staticmethod
  def sync_input(local_dir, remote_host, remote_dir, pattern):
    cmd_st = "ssh {0} mkdir -p {1}".format(remote_host, remote_dir)
    print cmd_st
    subprocess.Popen(cmd_st, shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    if pattern == None:
      pattern = ""

    cmd_st = "rsync --size-only -avtr {0}/{3} {1}:{2}".format(local_dir, remote_host, remote_dir, pattern)
    print cmd_st
    subprocess.Popen(cmd_st, shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()

  @staticmethod
  def sync_output(local_dir, remote_host, remote_dir, pattern):
    cmd_st = "mkdir -p {0}".format(local_dir)
    print cmd_st
    subprocess.Popen(cmd_st, shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()

    if pattern == None:
      pattern = ""
    cmd_st = "rsync --size-only -avtr {1}:{2}/{3} {0}".format(local_dir, remote_host, remote_dir, pattern)
    print cmd_st
    subprocess.Popen(cmd_st, shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()

#=====================================================================================

class LocalResource(Resource):
  
  default_config = { 'num_nodes': 1,
                     'job_concurrency': 1,
                     'num_deployments': 50,
                   }

  res_configs = { 'gpu': {'job_concurrency':  1, 'num_deployments': 8},
                  'dev': {'job_concurrency':  1, 'num_deployments': 1},
                }

  @staticmethod
  def get_paths():
    if 'Ubuntu' in platform.platform():
      oe_arch_dir = 'openeye/arch/Ubuntu-12.04-x64/oedocking/3.0.1'
    else:
      oe_arch_dir = 'centos/openeye/bin'
    path_dict = \
      {'charmm_bin':       '/damsl/projects/molecules/software/MD/Charmm/c36b1/exec/gnu/charmm',
       'namd_bin':        '/damsl/projects/molecules/software/MD/NAMD/NAMD/NAMD_2.9_Linux-x86_64-multicore/namd2',
       'gromacs_pdb2gmx': '/damsl/projects/molecules/software/MD/Gromacs/gromacs/gromacs-4.5.5/gromacs/bin/pdb2gmx',
       'gromacs_grompp':  '/damsl/projects/molecules/software/MD/Gromacs/gromacs/gromacs-4.5.5/gromacs/bin/grompp',
       'gromacs_mdrun':   '/damsl/projects/molecules/software/MD/Gromacs/gromacs/gromacs-4.5.5/gromacs/bin/mdrun',
       'omega2_bin':      '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/omega2',
       'fred_bin':        '/damsl/projects/molecules/software/DD/OpenEye/{0}/fred'.format(oe_arch_dir),
       'local_prefix':    Resource.local_prefix,
       'io_dir':          Resource.io_dir,
       'resource_prefix': '/damsl/projects/',
       'template_prefix': os.path.join('/home/{0}'.format(os.environ['USER']), '{0}/mddb/templates'),
       'template_link':   os.path.join('/home/{0}'.format(os.environ['USER']), '{0}/mddb/templates'),
       'data_dir':        '/damsl/mddb/data',
       'tleap':           '/damsl/projects/molecules/software/MD/Amber/amber12/bin/tleap',
       'sander':          '/damsl/projects/molecules/software/MD/Amber/amber12/bin/sander',
       'cuda_spfp':       '/damsl/projects/molecules/software/MD/Amber/amber12/bin/pmemd.cuda_SPFP',
       'AMBERHOME':       '/damsl/projects/molecules/software/MD/Amber/amber12/',
       'charmm':          '/damsl/projects/molecules/software/MD/Charmm/c36b1/exec/gnu/charmm',
       'charmm-medium':   '/damsl/projects/molecules/software/MD/Charmm/c36b1/exec/gnu/charmm',
       'namd':            '/damsl/projects/molecules/software/MD/NAMD/NAMD/namd2',
       'molcharge_bin':   '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/molcharge',
       'chimera':         '/damsl/projects/molecules/software/DD/Chimera/exec/bin/chimera',
       'sphgen':          '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/sphgen',
       'showbox':         '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/showbox',
       'grid':            '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/grid',
       'molcharge':       '/damsl/projects/molecules/software/DD/OpenEye/openeye/arch/redhat-RHEL5-x86/quacpac/1.6.3.1/bin/_molcharge-1.6.3.1',
       'omega2':          '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/omega2',
       'dock6':           '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/dock6',
       'sphere_selector': '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/sphere_selector',
       'apopdb2receptor': '/damsl/projects/molecules/software/DD/OpenEye/openeye/arch/redhat-RHEL6-x64/oedocking/3.0.1/bin/apopdb2receptor',
       'OE_LICENSE': '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/oe_license.txt'

      }

    return path_dict

  # overiding the remote deployment definition
  # The local deployment executes the generator directly in this method
  # Since the call is blocking, there is no need for probing for expected output files
  # like the remote deployment counterpart.
  def deploy_and_wait(self, session_dir, deployment_id, input_data, param_dict):
    deployment_dir            = self.get_deployment_name(deployment_id)
    abs_deployment_dir        = os.path.join(Resource.local_prefix, Resource.io_dir,  
                                             session_dir, deployment_dir)

    if not os.path.exists(abs_deployment_dir):
      os.makedirs(abs_deployment_dir)

    input_data_fn             = os.path.join(abs_deployment_dir, 'data.txt')
    with open(input_data_fn, 'w') as ofp:
      ofp.write(json.dumps(input_data))


    path_dict = LocalResource.get_paths()
    path_dict['template_prefix'] = path_dict['template_prefix'].format(self.user)
    path_dict['template_link']   = path_dict['template_link'].format(self.user)

    for d in input_data:
      #print d
      gen_name  = d.get('generator') or d.get('simulator') or d.get('docker')
      gen_class = Resource.generator_options[gen_name]
      gen_class.preprocess(d, path_dict)
      #print 'input_dict:', d
    #print 'input_data:', input_data

    self.execute_jobs_parallel(input_data, path_dict)
    conn = mddb_utils.get_dbconn(param_dict['dbname'],
                                 param_dict['dbuser'],
                                 param_dict['dbhost'],
                                 param_dict['dbpass'])


    self.load(conn, input_data, path_dict)
    conn.close()

    #self.cleanup(input_data, LocalResource.get_paths())

 

  @staticmethod
  def get_environ():
    d = { "PYTHONPATH": ["/scratch/01654/twoolf/projects/scripts/mddb/scheduler/:",
                         "/scratch/01654/twoolf/projects/scripts/mddb/simulation/:",
                         "/opt/apps/python/2.7.1/modules/lib/python:/opt/apps/python/2.7.1/lib:",
                         "/opt/apps/python/2.7.1/lib/python2.7/",
                        ],
          "PATH": ["/damsl/projects/molecules/software/MD/Charmm/c36b1/exec/gnu/",
                   "/opt/greenplum-db/./bin",
                   "/opt/greenplum-db/./ext/python/bin",
                   "/usr/local/bin",
                   "/bin",
                   "/usr/bin",
                   "/usr/local/sbin",
                   "/usr/sbin",
                   "/sbin"
                  ],
          "LD_LIBRARY_PATH": ["/lib64",
                              "/damsl/mddb/software/MD/Gromacs/centos/cpu-4.5.6/lib/",
                              "/opt/greenplum-db/./lib",
                              "/opt/greenplum-db/./ext/python/lib",
                              "/damsl/software/gearman/centos/lib",
                              "/damsl/software/boost/centos/lib",
                             ],
          "OE_LICENSE": ["/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/oe_license.txt"]
        }


    return d


   
#=====================================================================================
#
# A generic class for resources that use PBS-based systems.

class PBSResource(Resource):
  # Job submission method
  #   (1) Create an input data file, job script file, and resource configuration files
  #   (2) Copy these files to the remote resource
  #   (3) ssh into the remote resource and execute the appropriate job submission command
  def submit(self, session_dir, deployment_id, input_data):
    deployment_dir            = self.get_deployment_name(deployment_id)
    path_dict                 = self.__class__.get_paths()
    remote_prefix             = path_dict['resource_prefix']
    abs_deployment_dir        = os.path.join(Resource.local_prefix, Resource.io_dir,  
                                             session_dir, deployment_dir)
    abs_deployment_dir_remote = os.path.join(remote_prefix, Resource.io_dir, 
                                             session_dir, deployment_dir)
    input_data_fn             = os.path.join(abs_deployment_dir, 'data.txt')
    input_data_fn_remote      = os.path.join(abs_deployment_dir_remote, 'data.txt')
    jobscript_fn              = os.path.join(abs_deployment_dir, 'script.txt')
    jobscript_fn_remote       = os.path.join(abs_deployment_dir_remote, 'script.txt')
    res_configs_fn            = os.path.join(abs_deployment_dir, 'res.txt')
    res_configs_fn_remote     = os.path.join(abs_deployment_dir_remote, 'res.txt')

    with open(res_configs_fn, 'w') as ofp:
      ofp.write(json.dumps(self.res_configs))

    with open(input_data_fn, 'w') as ofp:
      ofp.write(json.dumps(input_data))

    with open(jobscript_fn, 'w') as ofp:
      ofp.write(self.compose_job_script(input_data_fn_remote, deployment_id, res_configs_fn_remote))

    self.__class__.sync_input_dirs(self.__class__.gateway_host, 
                                     [(os.path.join(Resource.io_dir, session_dir, deployment_dir), "*")], 
                                     Resource.local_prefix, 
                                     remote_prefix)
    cmd_st = "ssh {0} {1} {2}".format(self.gateway_host, self.submission_cmd, jobscript_fn_remote)
    print cmd_st
    subprocess.Popen(cmd_st,
                     shell=True, stdout=subprocess.PIPE,
                     stdin=subprocess.PIPE).stdout.read()

  def compose_job_script(input_data_fn_remote, deployment_id, res_configs_fn_remote):
    raise NotImplementedError( "Should have implemented this" )

#=====================================================================================

class LonestarResource(PBSResource):
  gateway_host = 'lonestar'
  submission_cmd = 'qsub'
  num_cores_per_node = 12

  default_config = {'num_nodes': 1,           # usually 1 (if not MPI)
                    'job_concurrency': 12,    # usually the number of cores per node (if not MPI)
                    'num_deployments': 10,    # max number of deployments: 50 (except for
                                              # the development queue, which is 2.
                    'time_limit': '24:00:00', # max time limit: 24 hours
                    'qname': 'normal',        # qname: normal, gpu, development
                    'sync_period': 3000,      # should be set to a minimum of 900 secs or
                                              # even longer if you want to be considerate. 
                   }

  res_configs = { 'mpi': {'num_nodes': 5, 'job_concurrency':  1, 'num_deployments': 10},
                  'gpu': {'job_concurrency':  2, # each Lonestar GPU node has 2 GPUs so 2 is the max here
                          'qname': 'gpu'},
                  'dev': {'qname': 'development', 'time_limit': '01:00:00', 
                          'sync_period': 900, # set it to the minimum requirement here for a good turn around time
                          'job_concurrency':6, 'num_deployments': 2},
                }


  # Compose a job script for Lonestar according to the resource configurations
  # and the specified input data file
  def compose_job_script(self, input_data_fn_remote, deployment_id, res_configs_fn_remote):
    qname = self.qname or LonestarResource.qname
    cmd = "/opt/apps/python/2.7.1/bin/python " +\
          "/scratch/01654/twoolf/projects/scripts/mddb/scheduler/resources.py " +\
          "--resource lonestar --mode execute " +\
          "--jobdata {0} --res_config_name {1} --res_configs_fn {2}"

    qsub_dict = { "num_cores": self.num_nodes * self.num_cores_per_node
                 ,"QUEUE": qname, "TIME": self.time_limit
                 ,"OUTFN": input_data_fn_remote + ".log"
                 ,"CMD": cmd.format(input_data_fn_remote, self.res_config_name, res_configs_fn_remote)
                }

    template = "#$ -S /bin/bash\n" +\
               "#$ -V\n" +\
               "#$ -j y\n" +\
               "#$ -o {OUTFN}\n" +\
               "#$ -pe 12way {num_cores}\n" +\
               "#$ -q {QUEUE}\n" +\
               "#$ -l h_rt={TIME}\n" +\
               "env " + Resource.print_environ(LonestarResource.get_environ()) + " {CMD}\n"

    return template.format(**qsub_dict)

  #def compose_ibrun_call(self, ibrun_dict):
  #  template = "ibrun -n {NUMCORES} -o {COREOFFSET} {EXECBIN} {INPARGS} > {OUTPUTLOG}"
  #  return template.format(**ibrun_dict)

  @staticmethod
  def get_environ():
    d = { "PYTHONPATH": ["/scratch/01654/twoolf/projects/scripts/mddb/scheduler/:",
                         "/scratch/01654/twoolf/projects/scripts/mddb/simulation/:",
                         "/opt/apps/python/2.7.1/modules/lib/python:/opt/apps/python/2.7.1/lib:",
                         "/opt/apps/python/2.7.1/lib/python2.7/",
                        ],
          "PATH": ["/opt/apps/python/2.7.1/bin",
                   "/opt/apps/tar-lustre/1.22/bin",
                   "/opt/apps/gzip-lustre/1.3.12/bin",
                   "/opt/apps/intel11_1/mvapich2/1.6/bin",
                   "/opt/apps/intel/11.1/bin/intel64",
                   "/opt/sge6.2/bin/lx24-amd64",
                   "/usr/lib64/qt-3.3/bin", 
                   "/usr/kerberos/bin",
                   "/usr/bin",
                   "/bin",
                   "/usr/sbin",
                   "/sbin",
                   "/opt/gsi-openssh-4.3/bin",
                   "/usr/X11R6/bin",
                   "/opt/ofed/bin",
                   "/opt/ofed/sbin",
                   "/usr/local/bin",
                   "/sge_common/default/pe_scripts",
                   ".",
                   "/work/01654/twoolf/Dock6/dock6/bin",
                  ],
          "LD_LIBRARY_PATH": ["/opt/apps/python/2.7.1/lib",
                              "/opt/apps/intel11_1/mvapich2/1.6/lib",
                              "/opt/apps/intel11_1/mvapich2/1.6/lib/shared",
                              "/opt/apps/intel/11.1/lib/intel64",
                              "/opt/gsi-openssh-4.3/lib",
                              "/opt/apps/gcc/4.4.5/lib64",
                             ],
          "AMBERHOME": ["/opt/apps/intel11_1/mvapich2_1_6/amber/12.0/"]
        }


    return d

  @staticmethod
  def get_paths():
    path_dict = { 'fred_bin':        '/work/01654/twoolf/oe/openeye/bin/fred'
                 ,'omega2_bin':      '/work/01654/twoolf/oe/openeye/bin/omega2'
                 ,'local_prefix':    Resource.local_prefix
                 ,'resource_prefix': '/scratch/01654/twoolf/projects/'
                 ,'io_dir':          Resource.io_dir
                 ,'charmm':          '/home1/01654/twoolf/c36a2_xlarge_dims_lonestar'
                 ,'charmm-medium':   '/work/01654/twoolf/charmm-medium'
                 ,'data_dir':        '/home1/01654/twoolf/data'
                 ,'namd':        '/home1/00288/tg455591/NAMD_2.8/NAMD_2.8_Linux-x86_64-MVAPICH-Intel-Lonestar/namd2'
                 ,'namd_mpi':    '/home1/00288/tg455591/NAMD_2.8/NAMD_2.8_Linux-x86_64-ibverbs-Lonestar/charmrun +p{nc} ++verbose ++mpiexec ++remote-shell /home1/00288/tg455591/NAMD_scripts/mpiexec ++runscript tacc_affinity /home1/00288/tg455591/NAMD_2.8/NAMD_2.8_Linux-x86_64-ibverbs-Lonestar/namd2'#.format(nc = LonestarResource.num_cores)
                 ,'template_prefix': '/scratch/01654/twoolf/projects/scripts/mddb/templates'
                 ,'template_link':   '/home1/01654/twoolf/templates' # don't forget to 'ln -s {template_prefix} {template_link}'. we use this when the template prefix is too long and not accepted by some programs written in fortran
                 ,'tleap':           '/opt/apps/intel11_1/mvapich2_1_6/amber/12.0/bin/tleap'
                 ,'sander':          '/opt/apps/intel11_1/mvapich2_1_6/amber/12.0/bin/sander'
                 ,'cuda_spfp':       '/opt/apps/intel11_1/mvapich2_1_6/amber/12.0/bin/pmemd.cuda_SPFP'
                 ,'AMBERHOME':       '/opt/apps/intel11_1/mvapich2_1_6/amber/12.0/'
                 ,'sphgen':       '/work/01654/twoolf/Dock6/dock6/bin/sphgen'
                 ,'showbox':       '/work/01654/twoolf/Dock6/dock6/bin/showbox'
                 ,'grid':       '/work/01654/twoolf/Dock6/dock6/bin/grid'
                 ,'dock6':       '/work/01654/twoolf/Dock6/dock6/bin/dock6'
                 ,'sphere_selector':       '/work/01654/twoolf/Dock6/dock6/bin/sphere_selector'
                 ,'chimera':       '/work/01654/twoolf/Chimera/exec/bin/chimera'
                 ,'molcharge':       '/work/01654/twoolf/oe/openeye/arch/redhat-RHEL5-x86/quacpac/1.6.3.1/molcharge'
                 ,'omega2':       '/work/01654/twoolf/oe/openeye/bin/omega2'
                }

    return path_dict
#=====================================================================================

class StampedeResource(PBSResource):
  gateway_host = 'stampede'
  num_cores_per_node = 16
  submission_cmd = 'sbatch'


  default_config = { 'num_nodes': 1, 
                     'job_concurrency': 16, 
                     'num_deployments': 10, 
                     'time_limit': '36:00:00',
                     'qname': 'normal',
                     'sync_period': 1800,
                   }

  res_configs = { 'mpi': {'num_nodes': 6, 'job_concurrency':  1, 'num_deployments': 10},
                  'gpu': {'job_concurrency':  1, 'qname': 'gpu'},
                  'dev': {'qname': 'development', 'time_limit': '04:00:00', 'sync_period': 900,
                          'job_concurrency': 1, 'num_deployments': 1
                         },
                }



  # The grammar is different from lonestar so we need a separate definition here
  def compose_job_script(self, input_data_fn_remote, deployment_id, res_configs_fn_remote):
    qname = self.qname or self.__class__.qname
    cmd = "/opt/apps/python/epd/7.3.2/bin/python " +\
          "/scratch/01654/twoolf/projects/scripts/mddb/scheduler/resources.py " +\
          "--resource stampede --mode execute " +\
          "--jobdata {0} --res_config_name {1} --res_configs_fn {2}"



    script_dict = { "num_cores": self.num_nodes * self.num_cores_per_node,
                    "JOBID": "{0:06x}".format(deployment_id),
                    "QUEUE": qname, "TIME": self.time_limit,
                    "OUTFN": input_data_fn_remote + ".log",
                    "CMD": cmd.format(input_data_fn_remote, self.res_config_name, res_configs_fn_remote),
                  }

    template = "#!/bin/bash\n" +\
               "#SBATCH -J {JOBID}         \n" +\
               "#SBATCH -o {OUTFN}         \n" +\
               "#SBATCH -n {num_cores}     \n" +\
               "#SBATCH -p {QUEUE}         \n" +\
               "#SBATCH -t {TIME}          \n" +\
               "env " + Resource.print_environ(StampedeResource.get_environ()) + " {CMD}\n"
    return template.format(**script_dict)

  #def compose_ibrun_call(self, ibrun_dict):
  #  template = "ibrun -n {NUMCORES} -o {COREOFFSET} {EXECBIN} {INPARGS} > {OUTPUTLOG}"
  #  return template.format(**ibrun_dict)

  @staticmethod
  def get_environ():
    d = { "PYTHONPATH": ["/scratch/01654/twoolf/projects/scripts/mddb/scheduler/:",
                         "/scratch/01654/twoolf/projects/scripts/mddb/simulation/:",
                         "/opt/apps/python/2.7.1/modules/lib/python:/opt/apps/python/2.7.1/lib:",
                         "/opt/apps/python/2.7.1/lib/python2.7/",
                        ],
          "PATH": ["/opt/apps/python/epd/7.3.2/bin/",
                   "/opt/apps/intel13/mvapich2/1.9/bin",
                   "/opt/apps/intel/13/vtune_amplifier_xe_2013/bin64",
                   "/opt/apps/intel/13/composer_xe_2013.2.146/mpirt/bin/intel64",
                   "/opt/intel/sep/bin",
                   "/usr/lib64/qt-3.3/bin",
                   "/usr/local/bin",
                   "/bin",
                   "/usr/bin",
                   "/usr/X11R6/bin",
                   "/opt/ofed/bin",
                   "/opt/ofed/sbin",
                   ".",
                   "/opt/apps/intel/13/composer_xe_2013.2.146/bin/intel64",
                   "/opt/apps/intel/13/composer_xe_2013.2.146/bin/intel64_mic",
                   "/opt/apps/intel/13/composer_xe_2013.2.146/debugger/gui/intel64",
                   "/work/01654/twoolf/oe/openeye/bin",
                  ],
          "LD_LIBRARY_PATH": ["/opt/apps/python/epd/7.3.2/lib",
                              "/opt/apps/cuda/5.0/lib64",
                              "/usr/lib64",
                              "/lib64",
                              "/opt/apps/intel13/mvapich2/1.9/lib",
                              "/opt/apps/intel13/mvapich2/1.9/lib/shared",
                              "/opt/apps/intel/13/composer_xe_2013.2.146/tbb/lib/intel64",
                              "/opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64",
                              "/opt/intel/mic/coi/host-linux-release/lib",
                              "/opt/intel/mic/myo/lib",
                              "/opt/apps/intel/13/composer_xe_2013.2.146/mpirt/lib/intel64",
                              "/opt/apps/intel/13/composer_xe_2013.2.146/ipp/../compiler/lib/intel64",
                              "/opt/apps/intel/13/composer_xe_2013.2.146/ipp/lib/intel64",
                              "/opt/apps/intel/13/composer_xe_2013.2.146/compiler/lib/intel64",
                              "/opt/apps/intel/13/composer_xe_2013.2.146/mkl/lib/intel64",
                              "/opt/apps/intel/13/composer_xe_2013.2.146/tbb/lib/intel64",
                             ],
          "AMBERHOME": ["/opt/apps/intel13/mvapich2_1_9/amber/12.0/"]
        }


    return d

  @staticmethod
  def get_paths():
    path_dict = { 'local_prefix':    Resource.local_prefix
                 ,'resource_prefix': '/scratch/01654/twoolf/projects/'
                 ,'io_dir':          Resource.io_dir
                 ,'charmm':          '/home1/01654/twoolf/c36a2_xlarge'
                 ,'charmm-medium':   '/home1/01654/twoolf/c36a2_xlarge'
                 ,'data_dir':        '/home1/01654/twoolf/data'
                 ,'namd':        '/home1/00288/tg455591/NAMD_2.8/NAMD_2.8_Linux-x86_64-MVAPICH-Intel-Lonestar/namd2'
                 ,'namd_mpi':    '/home1/00288/tg455591/NAMD_2.9_Source/Linux-x86_64-icc.ibverbs/charmrun +p{nc} ++verbose ++mpiexec ++remote-shell /home1/00288/tg455591/NAMD_scripts/mpiexec ++runscript tacc_affinity /home1/00288/tg455591/NAMD_2.9_Source/Linux-x86_64-icc.ibverbs/namd2'#.format(nc = LonestarResource.num_cores)

                 ,'template_prefix': '/scratch/01654/twoolf/projects/scripts/mddb/templates'
                 ,'template_link':   '/home1/01654/twoolf/templates'
                 ,'tleap':           '/opt/apps/intel13/mvapich2_1_9/amber/12.0/bin/tleap'
                 ,'sander':          '/opt/apps/intel13/mvapich2_1_9/amber/12.0/bin/sander'
                 ,'cuda_spfp':       '/opt/apps/intel13/mvapich2_1_9/amber/12.0/bin/pmemd.cuda_SPFP'
                 ,'AMBERHOME':       '/opt/apps/intel13/mvapich2_1_9/amber/12.0/'
                 ,'sphgen':       '/work/01654/twoolf/dock6/dock6/bin/sphgen'
                 ,'showbox':       '/work/01654/twoolf/dock6/dock6/bin/showbox'
                 ,'grid':       '/work/01654/twoolf/dock6/dock6/bin/grid'
                 ,'dock6':       '/work/01654/twoolf/dock6/dock6/bin/dock6'
                 ,'sphere_selector':       '/work/01654/twoolf/dock6/dock6/bin/sphere_selector'
                 ,'chimera':       '/work/01654/twoolf/Chimera/exec/bin/chimera'
                 #,'molcharge':       '/work/01654/twoolf/oe/openeye/bin/molcharge'
                 #,'omega2':       '/work/01654/twoolf/oe/openeye/arch/redhat-RHEL6-x64/omega/2.5.1.4'
 

                }

    return path_dict


#=====================================================================================
#*************************************************************************************
#*************************************************************************************
#*************************************************************************************
#*************************************************************************************
#*************************************************************************************

# Obsolete code below

#*************************************************************************************
#*************************************************************************************
#*************************************************************************************
#*************************************************************************************
#=====================================================================================
# A Cluster with a gateway node and a bunch or other nodes to directly run stuff on 
#
# A difference between a cluster resource and a PBS resource is that a cluster resource
# does not have a job scheduler running on the hear node. Compute jobs have to be
# deployed at a cluster compute node directly via ssh. That's why we need a list of
# compute nodes here.

class ClusterResource(Resource):
  gateway_host = None
  compute_nodes = []

  def submit(self, session_dir, deployment_dir, input_data):
    sync_list = []
    for input_dict in input_data:
      input_dict['session_dir'] = session_dir
      gen_name  = input_dict.get('generator') or input_dict.get('simulator') or input_dict.get('docker')
      gen_class = Resource.generator_options[gen_name]
      sync_list = sync_list + gen_class.get_sync_info(input_dict)

    sync_list = list(set(sync_list))

    path_dict            = self.__class__.get_paths()
    session_dir          = os.path.join(Resource.io_dir, session_dir)
    remote_prefix        = path_dict['resource_prefix']
    uid                  = 'job_' + str(uuid.uuid4())
    input_data_fn        = os.path.join(Resource.local_prefix,  session_dir, uid)
    input_data_fn_remote = os.path.join(remote_prefix, session_dir, uid)

    node_list            = self.__class__.compute_nodes
    num_nodes            = len(node_list)

    compute_node         = node_list[0]

    with open(input_data_fn, 'w') as ofp:
      ofp.write(json.dumps(input_data))


    self.__class__.sync_input_dirs(self.__class__.gateway_host, 
                                     sync_list + [(session_dir, "*")], 
                                     Resource.local_prefix, remote_prefix)
    cmd_st = self.__class__.compose_exec_cmd(compute_node, input_data_fn_remote)
    print cmd_st
    subprocess.Popen(cmd_st,
                     shell=True, stdout=subprocess.PIPE,
                     stdin=subprocess.PIPE)

 


#=====================================================================================
class DatascopeResource(ClusterResource):
  gateway_host = 'dsc'
  compute_nodes = ["dsp011", "dsp012", "dsp013", "dsp014", "dsp015", 
                   "dsp016", "dsp017", "dsp018", "dsp020"]

  @staticmethod
  def get_paths():
    d = { 'fred_bin':        '/home/yahmad/molecules/software/DD/OpenEye/openeye/arch/redhat-RHEL5-x64/oedocking/3.0.1/fred'
         ,'omega2_bin':      '/home/yahmad/molecules/software/DD/OpenEye/openeye/arch/redhat-RHEL5-x64/omega/2.4.6/omega2'
         ,'dock6_bin':       '/home/yahmad/molecules/software/DD/Dock/Dock/dock6/bin/dock6'
         ,'namd':       '/home/yahmad/molecules/software/MD/NAMD/NAMD/NAMD_2.9_Linux-x86_64-multicore/namd2'
         ,'charmm':      '/home/yahmad/molecules/software/MD/Charmm/c36b1/exec/gnu/charmm-xl'
         ,'charmm-medium':      '/home/yahmad/molecules/software/MD/Charmm/c36b1/exec/gnu/charmm-xl'
         ,'resource_prefix': '/home/yahmad'
         ,'io_dir':          Resource.io_dir
         ,'data_dir':        '/home/yahmad/molecules/data'
         ,'local_prefix':    Resource.local_prefix
         ,'template_prefix': '/home/yahmad/scripts/mddb/templates'
         ,'template_link':   '/home/yahmad/scripts/mddb/templates'
         #,'PYTHONPATH':      '/home/yahmad/scripts/mddb/simulation:/home/yahmad/nutanong/mddb/scheduler'
         ,'echo':            'echo'
         ,'tleap':           '/home/yahmad/MD/Amber/amber12/bin/tleap'
         ,'sander':          '/home/yahmad/MD/Amber/amber12/bin/sander'
         ,'cuda_spfp':       '/home/yahmad/MD/Amber/amber12/bin/pmemd.cuda_SPFP'
         ,'AMBERHOME':       '/home/yahmad/MD/Amber/amber12'
        }
    return d

  @staticmethod
  def get_environ():
    d = { 'AMBERHOME': ['/home/yahmad/MD/Amber/amber12'],
          'PATH': ['/usr/sbin:/sbin',
                   '/home/yahmad/cuda-5.0/bin',
                   '/home/yahmad/MD/Amber/amber12/bin',
                   '/usr/lib64/qt-3.3/bin',
                   '/usr/local/sbin',
                   '/usr/sbin',
                   '/sbin',
                   '/usr/local/bin',
                   '/bin',
                   '/usr/bin',
                   '/usr/X11R6/bin'],
          'LD_LIBRARY_PATH': ['/home/yahmad/ncarey/lib',
                              '/home/yahmad/cuda-5.0/lib64/lib',
                              '/home/yahmad/cuda-5.0/lib',
                              '/home/yahmad/cuda-5.0/lib64',
                              '/home/yahmad/MD/Amber/amber12/AmberTools/lib64',
                              '/home/yahmad/MD/Amber/amber12/lib64',
                              '/home/yahmad/MD/Amber/amber12/lib',
                              '/home/yahmad/MD/Amber/amber12/AmberTools/lib',
                              '/home/yahmad/MD/Amber/amber12/AmberTools/src/lib',
                              '/home/yahmad/MD/Amber/amber12/src/lib'],
          'PYTHONPATH': ['/home/yahmad/scripts/mddb/scheduler', 
                         '/home/yahmad/scripts/mddb/simulation']}
    return d

  @staticmethod
  def compose_exec_cmd(compute_node, input_data_fn_remote):
    env_st = Resource.print_environ(DatascopeResource.get_environ())
    cmd_st = "ssh {0} env {1} python2.6 /home/yahmad/scripts/mddb/scheduler/resources.py --resource datascope --mode execute --jobdata {2}".format(compute_node, env_st, input_data_fn_remote)

    return cmd_st




#=====================================================================================

class HHPCResource(PBSResource):
  gateway_host = 'hhpc'

  @staticmethod
  def get_paths():
    d = { 'fred_bin':        '/home/yahmad/molecules/software/DD/OpenEye/openeye/arch/redhat-RHEL5-x64/oedocking/3.0.1/fred'
         ,'omega2_bin':      '/home/yahmad/molecules/software/DD/OpenEye/openeye/arch/redhat-RHEL5-x64/omega/2.4.6/omega2'
         ,'dock6_bin':       '/home/yahmad/molecules/software/DD/Dock/Dock/dock6/bin/dock6'
         ,'namd2_bin':       '/home/yahmad/molecules/software/MD/NAMD/NAMD/NAMD_2.9_Linux-x86_64-multicore/namd2'
         ,'charmm':      '/home/yahmad/molecules/software/MD/Charmm/c36a2_xlarge_dims'
         ,'charmm-medium':      '/home/yahmad/molecules/software/MD/Charmm/c36a2_xlarge_dims'
         ,'resource_prefix': '/home/yahmad'
         ,'io_dir':          Resource.io_dir
         ,'data_dir':        '/home/yahmad/molecules/data'
         ,'local_prefix':    Resource.local_prefix
        }
    return d

  @staticmethod
  def compose_job_script(qname, input_data_fn_remote):
    qsub_dict = { "PE": "8way 96", "QUEUE": "batch", "TIME": "00:60:00"
                 ,"OUTFN": input_data_fn_remote + ".log"
                 ,"PYTHONPATH": '/home/yahmad/nutanong/mddb/simulation:/home/yahmad/nutanong/mddb/scheduler:/home/yahmad/software/site-packages'
                 ,"CMD": "python2.6 /home/yahmad/nutanong/mddb/scheduler/resources.py --resource hhpc --mode execute --jobdata {0}".format(input_data_fn_remote)}

    template = "#!/bin/bash\n" +\
               "#PBS -V\n" +\
               "#PBS -o {OUTFN}\n" +\
               "#PBS -e {OUTFN}.err\n" +\
               "#PBS -q {QUEUE}\n" +\
               "#PBS -l walltime={TIME}\n" +\
               "#PBS -l nodes=1\n" +\
               "cd /home/yahmad/nutanong/mddb\n" +\
               "env PYTHONPATH={PYTHONPATH} {CMD}\n"

    return template.format(**qsub_dict)
#*************************************************************************************
#*************************************************************************************
#*************************************************************************************
#*************************************************************************************
# Below is fine


#=====================================================================================
def test_simulators():
  rl = [LocalResource]



  for c in rl:
    print str(c),':'
    print '  Subclass:', issubclass(c, Resource)
    print '  Instance:', isinstance(c(), Resource)
    r = c()
    protein_seq = 'AAA'
    job_dict = \
      {'simulator':      'charmm'
       ,'trj_id':        1
       ,'protein_seq':   protein_seq
       ,'protein_id':    1
       ,'prefix':        os.path.join(r.path_prefix(),'data')
       ,'psf_file_name': None
       ,'sim_ff':        ''
       ,'phi_psi_array': []
      }
    d = r.deploy_and_wait(job_dict)
    print d

def test_dockers():
  rl = [LonestarResource]

  for c in rl:
    print str(c),':'
    print '  Subclass:', issubclass(c, Resource)
    print '  Instance:', isinstance(c(), Resource)
    r = c()

    #base_dir = '/scratch/01654/twoolf/projects/molecules/data/Odorant_GPCR'
    job_dict = \
      {'generator':           'fred'
       ,'job_id':             25
       ,'docker_id':          1
       ,'r_id':               1393
       ,'m_id':               1
       ,'c_id':               13
       ,'taut':               False
       ,'run_at':             'molecules/data/Odorant_GPCR/'
       ,'recepmod_dir':       'molecules/data/Odorant_GPCR/receptor_models/olfr1392/oe_fred_high_itasser_model1/'
       ,'compound_dir':       'molecules/data/Odorant_GPCR/compounds/'
       ,'recepmod_fred_fn':   'Olfactory_receptor_model_ITasser_1392_1_grid.oeb.gz'
       ,'smi_fn':             'ODL00000013_10024-56-3.smi'
      }
    print json.dumps(job_dict)

    d = r.deploy_and_wait(job_dict)
    print d

# /scratch/01654/twoolf/projects/molecules/data/Odorant_GPCR/receptor_models

def add_parser_options(option_parser, option_dict):
  for opt_name, default_val in option_dict.iteritems():
    try:
      default_val,comment = default_val
    except:
      comment = opt_name
      pass

    option_parser.add_option(
      "--"+opt_name,
      type=type(default_val), dest=opt_name,
      default=default_val,
      help=comment, metavar="#"+opt_name.upper())


def parse_param_dict(option_parser):
  options,args = option_parser.parse_args()
  return dict(vars(options).items())


#====================================================================
# Global Structure: Hostname -> Resource clas dict
# Add new resources here
resource_dict = { 'lonestar':  LonestarResource,
                  'stampede':  StampedeResource, 
                  'localhost': LocalResource,
                  'mddb-gpu':  LocalResource, 
                  'mddb2':     LocalResource, 
                  'hhpc':      HHPCResource,
                  'datascope': DatascopeResource,
                }

if __name__ == '__main__':
  param_dict = { "resource": "",
                 "mode":     "",
                 "jobdata":  "",
                 "res_configs_fn": "",
                 "res_config_name": "",
               }

  option_parser = OptionParser()
  add_parser_options(option_parser, param_dict)

  param_dict = parse_param_dict(option_parser)


  if param_dict['mode'] == 'execute':
    with open(param_dict['jobdata'], 'r') as ifp:
      data = json.loads(ifp.read())

    try:
      with open(param_dict['res_configs_fn'], 'r') as ifp:
        d = eval(ifp.read)
    except:
      d = {}
    res_class = resource_dict[param_dict['resource']]
    r         = res_class('', param_dict['res_config_name'], **d)


    res_paths = r.__class__.get_paths()
    print 'res_paths', res_paths
    r.execute_jobs_parallel(data, res_paths)
  else:
    test_dockers()












