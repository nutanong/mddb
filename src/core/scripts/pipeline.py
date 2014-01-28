
from abc import ABCMeta, abstractmethod, abstractproperty
import os
import time
import random
import math
import cStringIO
import itertools
import subprocess
import shutil
import time
from multiprocessing import Process

# Generic Pipeline Class to execute a sequence of shell commands and function calls

class Pipeline(object):
  __metaclass__ = ABCMeta

  # Application dependent and has to be implemented for each pipeline
  @abstractmethod
  def __init__(self, param_dict):
    raise NotImplementedError( "Should have implemented this" )

  # Application dependent and has to be implemented for each pipeline
  @abstractmethod
  def run_substage(self, param_dict):
    raise NotImplementedError( "Should have implemented this" )

  # A run wrapper with timeout
  def run(self):
    timeout = None
    try:
      timeout = self.timeout
    except:
      pass

    sleeptime = 5
    p = Process(target=self.run_)
    p.start()
    if not timeout:
      p.join()
      status = 'normal'
    else:
      total_sleep = 0

      while total_sleep < timeout:
        time.sleep(sleeptime)
        total_sleep = total_sleep + sleeptime
        if not p.is_alive():
          p.join()
          status = 'normal'
          break

      if p.is_alive():
        p.terminate()
        status = 'timeout'
      else:
        p.join()
        status = 'normal'

    return status

    

  # Actual Run method:
  #  (1) process template/input files by performing dict substitution using self.param_dict
  #  (2) rsync additional data (copying without substitutions)
  #  (3) run stages where each stage is a list of substages: [(command/operation, input, output)]
  #  (4) upon a substage failure, the execution restarts as the first substage of the stage
  def run_(self):


    self.param_dict['template_prefix'] = self.param_dict['template_prefix'].format(
                                           os.getenv('HOME'),
                                           self.param_dict.get('user'))

    output_prefix   = os.path.join(self.param_dict['resource_prefix'], 
                                   self.param_dict['io_dir'], 
                                   self.param_dict['session_dir'],
                                   self.param_dict['deployment_dir'],
                                   self.param_dict['run_dir'])

    self.param_dict['output_prefix'] = output_prefix

    if not os.path.exists(output_prefix):
      os.makedirs(output_prefix)

    for aux_dir in self.param_dict.get('aux_dirs') or []:
      cmd = 'rsync -av {0}/{1}/{2} {3}'.format(self.param_dict['template_prefix'], 
                                           self.param_dict['template_dir'], 
                                           aux_dir, output_prefix)
      print cmd
      subprocess.Popen(cmd, shell=True, 
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()

    if hasattr(self, 'aux_fn_list'):
      for fn in self.aux_fn_list:
        s = os.path.join(self.param_dict['template_prefix'], self.param_dict['template_dir'], fn)
        t = os.path.join(output_prefix, fn)

        work_dir = os.path.split(t)[0]
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)


        print s, '->', t
        shutil.copy2(s,t)

    for entry in  list(itertools.chain(self.inp_fn_list, *self.stage_list)):
      if isinstance(entry, dict):
        subdir   = entry.get('d') or ""
        fns = [entry.get('i') or entry.get('si'), entry.get('a')]
      else:
        subdir   = ""
        fns = [entry]

      for fn in fns:
        if fn != None:
          with open(os.path.join(self.param_dict['template_prefix'],
                                 self.param_dict['template_dir'], subdir, fn), 'r') as ifp:
            template = ifp.read()

            fn = os.path.join(output_prefix, subdir, fn)
            work_dir = os.path.split(fn)[0]
            print 'work_dir:', work_dir
            if not os.path.exists(work_dir):
              os.makedirs(work_dir)


            with open(fn, 'w') as ofp:
              ofp.write(template.format(**self.param_dict))

    for stage in self.stage_list:
      move_on = False
      while not move_on:
        for substage in stage:
          out_fns = substage.get('o')
          in_fn   = substage.get('i') or substage.get('si')
          if out_fns == None:
            try:
              out_fns = [os.path.splitext(in_fn)[0] + self.default_out_ext]
              substage['o'] = out_fns
              cmd = 'rm ' + ' '.join(out_fns)
              print cmd
              subprocess.Popen(cmd, shell=True, cwd = output_prefix,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()



            except:
              substage['o'] = None
         
          success = self.run_substage(output_prefix, substage)
          #success = True
          if success:
            move_on = True
          else:
            move_on = False
            break
    print 'pipleline finished'

