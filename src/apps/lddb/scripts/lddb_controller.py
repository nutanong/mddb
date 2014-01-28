from optparse import OptionParser

import re
import StringIO

import datetime, time

import xml.etree.ElementTree as ET
import subprocess
import os
import sys
import random
import socket
import pprint
import mddb_utils
import json
import itertools
import inspect
import glob
from mddb_cmd import *

import parser as mddb_parser
import cStringIO


import math

from multiprocessing import Process

import controller
import psycopg2  
import psycopg2.extras
import resources

class LDDBController(controller.Controller):

# Uncomment this method if we need to add more param dict keys
#
#  def __init__(self):
#    self.param_dict['newkey'] = ''
#    self.add_parser_options(self.param_dict)
#    self.param_dict = self.parse_param_dict()
#    sys.stdout.flush()
#

  def run(self):
    mode_dict = {'initdb': self.initdb,
                 'work': self.work,
                 'load': self.load,
                 'loop': self.control_loop,
                 'load_trj': self.load_trj,
                }

    mode = self.param_dict['mode']
    mode_dict[mode]()

  # Initilize database:
  #  (1) quit worker processes; 
  #  (2) load database schema files;
  #  (3) load stored procedure files
  # We need to do this only once or whenever we want to wipe the data.
  # Example: python scheduler/lddb_controller.py --dbname lddb_odorant2 --mode initdb
  def initdb(self):
    comment = "mode: initdb\nRequired Param: --dbname"
    print comment
    schema_files = ['job_control.sql', 'lddb_schema.sql']
    self.quit_everything()
    self.init_database(schema_files)
    self.load_stored_procedures(['lddb_utils.sql'])
    
    self.load_compound_states()
    self.load_receptor_frames()

  # Work: Main control loop
  #   (1) reload stored procedures (in case some of them changed after initdb);
  #   (2) quit worker processes;
  #   (3) reinit the job control data
  #   (4) setup worker processes by specifying a list of dictionaries where
  #       each dictionary is associated with a particular resource and
  #       contains: the host on which data generation will take place, 
  #       average number of sequential jobs, a predefined set of configuration parameters
  #   (5) read unfinished jobs from the database and put them in the job queue for
  #       the worker processes to consume
  #   (6) (not implemented yet) start a control loop which keeps queueing more jobs
  #       according to some statistical analysis results coming back from previous jobs
  # Example: python scheduler/lddb_controller.py --dbname lddb_odorant2 --mode work
  def work(self):
    comment = "mode: initdb\nRequired Param: --dbname"
    print comment

    self.load_stored_procedures(['lddb_utils.sql'])
    self.quit_everything()
    self.load_stored_procedures(['job_control.sql'])
    l = []
    l.append({'hostname':'localhost', 'avg_seq_jobs': '1', 'res_config_name': 'def'})
    #l.append({'hostname':'stampede',  'avg_seq_jobs': '1', 'res_config_name': 'dev'})
    #l.append({'hostname':'lonestar',  'avg_seq_jobs': '1', 'res_config_name': 'mpi'})

    self.setup_workers(l)
    conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                 self.param_dict['dbuser'],
                                 self.param_dict['dbhost'],
                                 self.param_dict['dbpass'])

    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    cur.execute("select * from level1_unfinishedjobsview order by job_id;")
    #cur.execute("select * from level1_unfinishedjobsview order by job_id;")
    #cur.execute("select 'membrane' as generator,* from unprocessedreceptormodelsview order by receptor_id,model_id;")

   

    run_params = cur.fetchall()
    #run_params = run_params[:5]

    for d in run_params:
      d['dbname']    = self.param_dict['dbname']
      data = json.dumps(d)
      job_id = d.get('job_id') or 'NULL'
      cur.execute("select jobqueue_insert({0},'{1}')".format(job_id, data))

    cur.close()
    conn.commit()
    conn.close()

    # start a screen running the controller in the loop mode
    # cmd = python scheduler/lddb_controller --dbname {dbname} --mode loop

  def control_loop(self):
    while True:
      # run some analysis, e.g., a stored procedure inside the database
      # put more jobs in the job queue
      time.sleep(10)

  @staticmethod
  def get_recpmodel_dir(d):
    bname = os.path.splitext(d['pdb_fn'])[0]

    prefix = '/damsl/projects/molecules/data/Odorant_GPCR/receptor_models'
    n = os.path.join(prefix,
                     'olfr{0}'.format(d['receptor_id']),
                     d['docker'],
                     '{0}_{1}'.format(d['model_id'], bname),
                     'bilayer_frames')
    return n

  # If the worker process was killed before loading the data
  # we can use this mode to load the data in the run directory
  # by specifying 'jobdata' which is the path to the job dictionary
  # in the run directory
  def load(self):
    with open(self.param_dict['jobdata'], 'r') as ifp:
      data = json.loads(ifp.read())

    conn = psycopg2.connect(dbname=self.param_dict['dbname'])
    resources.Resource.load(conn, data, resources.LocalResource.get_paths())
    conn.close()

  def load_trj(self):
    jobdata = self.param_dict['jobdata']
    job_dir = os.path.split(jobdata)[0]
    print 'job_dir:', job_dir
    with open(jobdata, 'r') as ifp:
      d = json.loads(ifp.read())

    conn = psycopg2.connect(dbname=self.param_dict['dbname'])
    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    cur.execute('select receptormodel_id,pdb_fn from ReceptorModels where receptor_id = {receptor_id} and model_id = {model_id}'.format(**d))
    d.update(cur.fetchone())
    d['docker'] = 'dock6'
    target_dir = LDDBController.get_recpmodel_dir(d)
    frames_dir = os.path.join(target_dir, 'frames')
    if not os.path.exists(frames_dir):
      os.makedirs(frames_dir)
    
    print 'copying pdb and dcd from {0} to {1}'.format(job_dir, target_dir)
    cmd_st = 'cp step5_assembly.pdb step5_assembly.xplor.psf namd/step7.1_production.dcd {0}'.format(target_dir)
    #subprocess.call(cmd_st, shell = True, cwd = job_dir)

    cmd_st = '/damsl/projects/molecules/software/tools/VMD/bin/vmd -dispdev text -e /damsl/projects/molecules/data/Odorant_GPCR/receptor_models/extract_frames.tcl'
    #subprocess.call(cmd_st, shell = True, cwd = frames_dir)

    output = cStringIO.StringIO()
    output.write("{0},T,T,{1}\n".format(d['receptormodel_id'], frames_dir))
    output.seek(0)


    cur.execute("select count(*) from Trajectories")
    print cur.fetchone()['count']

    cur.copy_from(output, 'Trajectories', sep = ',', columns = ('receptormodel_id', 'expsolv', 'bilayer', 'trj_dir'))
    cur.execute("select trj_id from Trajectories where trj_dir = '{0}'".format(frames_dir))
    trj_id = cur.fetchone()['trj_id']
    print trj_id
    cur.execute("select count(*) from Trajectories")
    print cur.fetchone()['count']


    output.close()
    output = cStringIO.StringIO()
    for frame_fn in os.listdir(frames_dir):
      t = os.path.splitext(frame_fn)[0]
      t = int(t)
      output.write("{0},{1},{2}\n".format(trj_id, t, frame_fn))
    output.seek(0)
    cur.copy_from(output, 'Frames', sep = ',', columns = ('trj_id', 't', 'pdb_fn'))
    output.close()
    cur.close()
    conn.commit()

  def load_compounds(self):
    smi_list_fn = '/damsl/projects/molecules/data/Odorant_GPCR/smi_list.txt'
    with open(smi_list_fn, 'r') as ifp:
      smi_list = ifp.readlines()


    output = StringIO.StringIO()
    dict   = {}

    for f in smi_list:
      dir,f = os.path.split(f)
      m = re.search('ODL(.+?)_', f)
      c_id = int(m.group(1))
      t = 'tautomers' in f
      if not t:
        with open(os.path.join('/damsl/projects',dir, f.rstrip()), 'r') as smi_fp:
          smi_st = smi_fp.read().split()[0]

        output.write("{0},{1},{2},{3},{4}\n".format(c_id,t, f.rstrip(), smi_st, dir))

    output.seek(0)
    conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                 self.param_dict['dbuser'],
                                 self.param_dict['dbhost'],
                                 self.param_dict['dbpass'])
    cur = conn.cursor()
    cur.copy_from(output, 'InitialCompounds', sep = ',',
                  columns = ('compound_id', 'tautomers', 'smi_fn', 'smiles', 'compound_dir'))

    cur.execute('select count(*) from InitialCompounds')
    print cur.fetchone()[0], 'compounds added'
    cur.close
    conn.commit()
    conn.close()
    output.close()  

  def load_compound_states(self):
    self.load_compounds()

    prefix = '/damsl/projects/'
    taut_cmd = '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/tautomers'

    conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                 self.param_dict['dbuser'],
                                 self.param_dict['dbhost'],
                                 self.param_dict['dbpass'])

    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    cur.execute('select compound_dir,smi_fn,compound_id from initialcompounds order by compound_id')
    compound_dicts = cur.fetchall()
    #cur.execute('truncate table CompoundStates')

    for d in compound_dicts:
      output = StringIO.StringIO()
      cmd = taut_cmd + ' ' + os.path.join(prefix,d['compound_dir'], d['smi_fn'])
      taut_fn = '/tmp/tautout{0}.txt'.format(d['compound_id'])
      with open(taut_fn, 'w') as ofp:
        subprocess.call(cmd, shell=True, stdout=ofp,
                         cwd=os.environ['HOME']
                       )
      with open(taut_fn, 'r') as ifp:
        l = ifp.readlines()
        print l
        for state_id,compound_state_smi in enumerate(l):
          compound_state_smi = compound_state_smi.split()[0]
          line = '{0}\t{1}\t{2}\n'.format(d['compound_id'], state_id, compound_state_smi)
          print line
          output.write(line)

      output.seek(0)
      cur.copy_from(output, 'CompoundStates', columns = ('compound_id', 'state_id', 'smiles'))
      conn.commit()

    conn.close()


#
## ls -1 olfr*/*dock6*/*.bmp olfr*/*fred_high*/gpcrm_*_olfr????_grid.oeb.gz  olfr*/*itasser*/*.bmp
#
#
  def load_receptor_frames(self):
    data_dir = '/damsl/projects/molecules/data/Odorant_GPCR/'
    frames_fn         = os.path.join(data_dir, 'frames.txt')
    trajectories_fn   = os.path.join(data_dir, 'trajectories.txt')
    receptormodels_fn = os.path.join(data_dir, 'receptormodels.txt')
    conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                 self.param_dict['dbuser'],
                                 self.param_dict['dbhost'],
                                 self.param_dict['dbpass'])
    cur = conn.cursor()
    with open(receptormodels_fn, 'r') as ifp:
      cur.copy_from(ifp, 'ReceptorModels',
                    columns = ('receptor_id', 'model_id', 'pdb_fn'))
      cur.execute('select count(*) from ReceptorModels')
      print cur.fetchone()[0], 'receptor models added'

    
    with open(trajectories_fn, 'r') as ifp:
      cur.copy_from(ifp, 'Trajectories',
                    columns = ('receptormodel_id', 'expsolv', 'bilayer', 'trj_dir'))
      cur.execute('select count(*) from Trajectories')
      print cur.fetchone()[0], 'trajectories added'

    with open(frames_fn, 'r') as ifp:
      cur.copy_from(ifp, 'Frames',
                    columns = ('trj_id', 't', 'fred_fn', 'dock6_fn', 'pdb_fn'))
      cur.execute('select count(*) from Frames')
      print cur.fetchone()[0], 'frames added'



    cur.close
    conn.commit()
    conn.close()
  
if __name__ == '__main__':
  if LDDBController.check_user():
    controller = LDDBController()
    controller.run() 
  
  
