from optparse import OptionParser

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
import controller

import ast
import math

from multiprocessing import Process
import psycopg2
import psycopg2.extras
import chunkIt


try:
  mwd = os.environ['MDDBWORKDIR']
except:
  mwd = os.environ['HOME'] + "/mddb"

ubuntu_pythonpaths = "/damsl/software/python/ubuntu/lib/python2.7/site-packages/MDAnalysis-0.7.6-py2.7-linux-x86_64.egg/:/damsl/software/python/ubuntu/lib/python2.7/site-packages/:/usr/include/:/damsl/projects/molecules/software/tools/mddb/:{mddb_work_dir}/simulation/:{mddb_work_dir}/mddb/scheduler/".format(mddb_work_dir = mwd)

centos_pythonpaths = "/damsl/software/python/centos/lib/python2.6/site-packages/:/usr/include/:{mddb_work_dir}/scheduler/:{mddb_work_dir}/simulation/:/damsl/projects/molecules/software/tools/mddb/:/opt/greenplum-db/./lib/python".format(mddb_work_dir = mwd)


config_dir = '/damsl/mddb/data/mddb_configs'

def getStartingPoints(trj_id, trj_len, num_points):
  timestamps = []
  dict_list  = []
  while len(timestamps) < num_points:
    t = random.randint(0, trj_len - 1)
    if not t in timestamps:
      timestamps.append(t)
      dict_list.append({'trj_id': trj_id, 't': t})

  return dict_list
  


class MDDBController(controller.Controller):
  mddb_param_dict = {
    'gmw':           'scheduler/worker.py'.format(mddb_work_dir = mwd)
    ,'timestep_pico': 0.001
    ,'trj_save_freq': 100
    ,'jobfile':       ''
    ,'csv_prefix':   '/damsl/projects/molecules/data/bpti_csv/amb_amber'.format(mddb_work_dir = mwd)
    ,'protein_name': 'bpti'
    ,'config_id': 1
    ,'trj_list_fn':  'bpti-all_times.csv'
    ,'trj_dir': '/damsl/mddb/data/deshaw/MD0_DATA/Public/outgoing/science2010/DESRES-Trajectory-bpti-all/bpti-all'
    ,'pdb_dir': '/damsl/mddb/data/deshaw/MD0_DATA/Public/outgoing/science2010/DESRES-Trajectory-bpti-all/pdb_states'
    ,'pdb_fn': 'state1.pdb'
    ,'user_id': '1'
    ,'policy':  'long running'
  }
  
  schema_files = ['controller.sql', 'schema.sql']
  sp_files     = ['utils.sql', 'maintenance.sql', 'decision.sql', 'loading.sql', 'control.sql']

  @staticmethod
  def get_prev_timestamp(trj_id, t):
    t = t-1
    if t == -1:
      trj_id = trj_id -1
      t = 999

    return trj_id,t


  @staticmethod
  def get_prev_timestamps(trj_id, t, n):
    l = []
    
    for i in range(0, n):
      trj_id,t = MDDBController.get_prev_timestamp(trj_id, t)
      if trj_id < 1:
        break
      l.append((trj_id, t))

    return l

  def __init__(self):
    self.param_dict = dict(self.param_dict.items() + self.mddb_param_dict.items())
    self.add_parser_options(self.param_dict)
    self.param_dict = self.parse_param_dict()

    if self.param_dict['mode'] == 'initdb':
      self.quit_everything()
      self.init_database(self.schema_files)
      self.load_stored_procedures(self.sp_files)
      self.load_configs()
      sys.exit(0)

    if self.param_dict['mode'] == 'initss':
      mddb_utils.subspace_init(self.param_dict)
      sys.exit(0)

    if self.param_dict['mode'] == 'trjs':
      conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                   self.param_dict['dbuser'],
                                   self.param_dict['dbhost'],
                                   self.param_dict['dbpass'])

      mddb_utils.add_trajectories(conn, self.param_dict)
      conn.close()
      sys.exit(0)




    if self.param_dict['mode'] == 'loop':
      self.sample_loop()
      sys.exit(0)
  
    if self.param_dict['mode'] == 'setup':
      self.quit_everything()
      #self.clear_jobqueue()
      #self.init_database(self.schema_files)
      #self.load_stored_procedures(self.sp_files)
      #self.load_configs()
      l = []
      l.append({'hostname':'stampede', 'sequentialism': '1', 'res_config_name': 'gpu'})

      self.setup_workers(l)
      #self.start_sampleloop()
      sys.exit(0)

    if self.param_dict['mode'] == 'qjobs':
      #trj_id_list = ['137']
      #trj_len = 1000
      #num_points = 1 



      #l = map(lambda x: getStartingPoints(x, trj_len, num_points), trj_id_list)

      #l = list(itertools.chain(*l))

      l = []
      timestamps = [(3294, 555, 50), (3264, 5, 50), (3264, 929, 100)]
      for ts in timestamps:
        l = l + MDDBController.get_prev_timestamps(*ts)

      l = map(lambda (trj_id,t): {'trj_id': "{0:03d}".format(trj_id), 't': t}, l)

      l = l[99:]

      default_params = {'nstep_simulation': 50000000,
                        'trj_save_freq': 50000,
                        'source': 'deshaw',
                        'generator': 'amber',
                        'template_dir': 'amber_min',
                        'dbname': self.param_dict['dbname'],
                        'dbuser': self.param_dict['dbuser'],
                        'dbhost': self.param_dict['dbhost'],
                        'dbpass': self.param_dict['dbpass'],
                       }


      conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                   self.param_dict['dbuser'],
                                   self.param_dict['dbhost'],
                                   self.param_dict['dbpass'])


      cur = conn.cursor()
      cur.execute('truncate table jobqueue')
      for d in l:
        d.update(default_params)
        print d
        data = 'ARRAY' + str([json.dumps(d)])
        cur.execute('select jobqueue_insert({0})'.format(data))


      cur.close()
      conn.commit()
      conn.close()      


    if self.param_dict['mode'] == 'cmd':
      mc = mddb_cmd()
      mc.init_conn(self.param_dict)
      mc.cmdloop()
      sys.exit(0)
 
  def clear_jobqueue(self):
    conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                 self.param_dict['dbuser'],
                                 self.param_dict['dbhost'],
                                 self.param_dict['dbpass'])

    cur = conn.cursor()

    cur.execute('truncate table JobQueue')
    cur.close()
    conn.commit()
    conn.close()


  def quit_everything(self):
    mddb_utils.run_cmd('screen -S {0}_{1}_sample_loop -X quit'.format(self.param_dict['dbname'],
                                                                      self.param_dict['user']))
    super(MDDBController, self).quit_everything()

  def load_configs(self):
    conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                 self.param_dict['dbuser'],
                                 self.param_dict['dbhost'],
                                 self.param_dict['dbpass'])
  
    c_id_all = mddb_utils.execute_query(conn, "select insert_new_config('charmm', 'all')")
    self.load_config(c_id_all, '/damsl/mddb/data/amino_configs/config_charmm')
    c_id_ua = mddb_utils.execute_query(conn, "select insert_new_config('charmm', 'ua')")
    self.load_config(c_id_ua,  '/damsl/mddb/data/amino_configs/config_ua_charmm')
  
    #mddb_utils.load_machine_info(conn, ['mddb','mddb2','qp2','qp3','qp4','qp5','qp6'])
    conn.commit()
    conn.close() 
  
  def read_config_file(config_file):
    print "Parsing configuration file ", config_file
  
    tree = ET.parse(config_file)
    root = tree.getroot()
  
    for param_name in param_dict.iterkeys():
      node = root.find(param_name)
      if node != None:
        param_dict[param_name] = node.text
        print "Setting '{0}' to '{1}'".format(param_name, node.text)
        
    for k,v in param_dict.items():
      if v == '':
        print "{0}: Error '{1}' is NOT defined.".format(config_file, k)
        exit(1)
  
  def load_config(self, c_id, prefix):
    mddb_utils.run_sql("select load_config({0},'{1}')".format(
                  c_id, prefix), self.param_dict['dbname']);


  def sample_loop(self):
    while True:
      #try:
        conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                     self.param_dict['dbuser'],
                                     self.param_dict['dbhost'],
                                     self.param_dict['dbpass'])
    
    
        res = mddb_utils.execute_query(conn,
          "select param_value from ControlParams where param_name = 'sampleloop'")
    
        print datetime.datetime.now()
        try:
          nt = mddb_utils.execute_query(conn, "select resample_blocking(E'{gmdhost}:{gmdport}')".format(**self.param_dict))
        except Exception as e:
          m = "Error {0}".format(str(e))
          print m
          try:
            mddb_utils.execute_query(conn, "insert into ErrorMessages values (statement_timestamp(),{0},{1})".format("sample_loop", m.encode("utf-8")))
          except:
            pass
          nt = 0
          time.sleep(2)
          continue
  
        if int(nt) > 0:
          time.sleep(0.05)
        else:
          time.sleep(2)
  
        print res
        if res == 'quitting':
          mddb_utils.execute_query(conn,
            "select controlparams_upsert('sampleloop', 'notrunning')")
          mddb_utils.execute_query(conn,
            "update ExpJobs set active = False;")
    
          conn.close()
          break
    
        cur = conn.cursor()
        cur.execute('select S.subspace_id from subspaces S, expjobs E where S.subspace_id = E.subspace_id group by S.subspace_id having bool_or(active)')

        num_active_subspaces = cur.rowcount
        
        if num_active_subspaces == 0:
          cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
          cur.execute('select * from mdq_pop()')
          if cur.rowcount > 0:
            ss_dict = cur.fetchone()

            if any(ss_dict.values()):
              ss_dict['protein_seq'] = "".join(ss_dict['expanded_terms'])
              mddb_utils.subspace_init(dict(self.param_dict.items() + ss_dict.items()))

          cur.close()
          conn.commit()
          conn.close()
  
  def start_sampleloop(self):
    conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                 self.param_dict['dbuser'],
                                 self.param_dict['dbhost'],
                                 self.param_dict['dbpass'])
 
    prefix = os.path.join(os.getenv('HOME'),
                          self.param_dict['user'],
                          controller.Controller.mwd)


    pp = controller.Controller.centos_pythonpaths.format(prefix)
 
 
    env = "PYTHONPATH={0} GPHOME=/opt/greenplum-db/. PGHOST=localhost4 PGOPTIONS=--client-min-messages=error ".format(pp)
    env = env + "PATH={PATH} LD_LIBRARY_PATH={LD_LIBRARY_PATH} ".format(**os.environ)

    self.param_dict['env'] = env

    sample_loop_cmd = ("""screen -S {dbname}_{user}_sample_loop -d -m """+\
                       """env {env}""" +\
                       """/opt/greenplum-db/ext/python/bin/python """+ os.path.abspath( __file__ ) +\
                       """ --mode loop """+\
                       """ --dbname  {dbname} """+\
                       """ --dbpass  {dbpass} """+\
                       """ 2>&1 | tee /tmp/{dbname}_{user}_sampleloop.log >&1 """).format(**self.param_dict);
  
  
  
    mddb_utils.run_cmd(sample_loop_cmd)
  
    mddb_utils.execute_query(conn,
      "select controlparams_upsert('sampleloop', 'running')")
  
    conn.commit()
    conn.close()
  
  def queue_jobs(self, jobfile):
    conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                 self.param_dict['dbuser'],
                                 self.param_dict['dbhost'],
                                 self.param_dict['dbpass'])

    with open(jobfile, 'r') as ifp:
      content = ifp.read()
      print content
      job_obj = ast.literal_eval(content)  #json.loads(content, object_hook=mddb_utils.ascii_encode_dict)


      print job_obj
      job_dicts = controller.Controller.get_job_dicts(**job_obj)

    for d in job_dicts:
      print 'queue job: ', d

      mddb_utils.mdq_push(conn,**d)
      #mddb_utils.mdq_push(conn, 0, psp_terms, num_expjobs, num_iterations, num_steps_sim, u_id, policy,
      #                    simulator, config_id)

    return

if __name__ == '__main__':
  if MDDBController.check_user():
    controller = MDDBController()
