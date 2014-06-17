from optparse import OptionParser

import re
import StringIO

import datetime, time

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
import resources

from multiprocessing import Process
  
# Generic controller class
class Controller(object): 
  sharedacc = 'webapp'
  centos_pythonpaths = "/damsl/projects/molecules/software/SP/Modeller/modeller-9.11/modlib/:/damsl/projects/molecules/software/SP/Modeller/modeller-9.11/lib/x86_64-intel8/python2.5/:/damsl/software/python/centos/lib/python2.6/site-packages/:/usr/include/:{0}/scheduler/:{0}/simulation/:/damsl/projects/molecules/software/tools/mddb/:/opt/greenplum-db/./lib/python"

  ubuntu_pythonpaths = "/damsl/software/python/ubuntu/lib/python2.7/site-packages/MDAnalysis-0.7.6-py2.7-linux-x86_64.egg/:/damsl/software/python/ubuntu/lib/python2.7/site-packages/:/usr/include/:/damsl/projects/molecules/software/tools/mddb/:{0}/simulation/:{0}/mddb/scheduler/"

  mwd = 'mddb' 
  
  option_parser = OptionParser()

  # default parameter dictionary for command line arguments
  param_dict = {
    'mwd'         : mwd,
    'mode'        : 'cmd',
    'dbdir'       : 'database',
    'dbname'      : ('testdb',
                     'name of the database'),
    'dbhost'      : os.getenv('HOSTNAME') or socket.gethostname(),
    'dbuser'      : os.getenv('USER'),
    'gmdhost'     : os.getenv('HOSTNAME') or socket.gethostname(),
    'dbpass'      : '{0}/.mypass'.format(os.getenv('HOME')),
    'prefix'      : '/damsl/mddb/data',
    'gmd'         : ('/damsl/software/gearman/centos/sbin/gearmand',
                     'gearmand executable'),
    'gmw'         : ('scheduler/worker.py',
                     'gearman worker script'),
    'madpack'     : '/usr/local/madlib/bin/madpack',
    'jobfile'     : '',
    'user'        : '',
    'jobdata'     : '',
  }
 
  # generic __init__ method which parses command line parameters
  def __init__(self, additional_dict=None):
    if additional_dict != None:
      self.param_dict = dict(self.param_dict.items() + self.additional_dict.items())
    self.add_parser_options(self.param_dict)
    self.param_dict = self.parse_param_dict()
    sys.stdout.flush()


  @staticmethod
  def get_job_dicts(default_dict, range_dicts):
    if range_dicts == []:
      return [default_dict]

    res = []
    for r in range_dicts:
      list_of_list_of_pairs = []
      for k,v in r.items():
        pairs = []
        for a in v:
          pairs.append((k,a))
        list_of_list_of_pairs.append(pairs)

      new_res = map(lambda x: dict(default_dict.items() + dict(x).items()),
                    list(itertools.product(*list_of_list_of_pairs)))
      res = res + new_res
    return res   

  #iterates through the default param dict and add each entry into
  #command line argument options. For example, if you want to do the following
  # self.option_parser = OptionParser()
  # self.option_parser.add_option("--file", type=str, dest="file", default="",
  #                               help="write report to FILE", metavar="FILE")
  # self.option_parser.add_option("--user", type=str, dest="user", default="guest",
  #                                help="specify user name", metavar="USER")
  #you would do the following instead
  # self.add_parser_options({'file': ("","write report to FILE"), 
  #                          'user': ("guest","specifying user name")})
  # The default parameter dictionary is defined as self.default_dict above

  def add_parser_options(self,option_dict):
    for opt_name, default_val in option_dict.iteritems():
      try:
        default_val,comment = default_val
      except:
        comment = opt_name
        pass

      self.option_parser.add_option(
        "--"+opt_name, 
        type=type(default_val), dest=opt_name,
        default=default_val,
        help=comment, metavar="#"+opt_name.upper())

  #does the actual parsing of parameters
  def parse_param_dict(self):
    options,args = self.option_parser.parse_args()
    return dict(vars(options).items())
    
 

  # database init function:
  #   (1) try to drop the current database and to create a new one 
  #   (2) load schema files specified in 'init_files
  # call this method when you want to create a new database or wipe an existing one
  def init_database(self, init_files):
    dbname = self.param_dict['dbname']
    dbdir = os.path.join(os.getenv('HOME'),
                         self.param_dict['user'],
                         self.param_dict['mwd'],
                         self.param_dict['dbdir'])
    p = mddb_utils.run_cmd("dropdb {0}".format(dbname))
    p = mddb_utils.run_cmd("createdb {0}".format(dbname))

    madpack    = '/usr/local/madlib/bin/madpack'
    #mddb_utils.run_sql("alter role {0} with password '12345'".format(self.param_dict['dbuser']), dbname)
    #mddb_utils.run_cmd("{0} -p greenplum -c {1}/12345@{2}/{3} install".format(
    #          madpack, self.param_dict['dbuser'], 'localhost4', dbname))
  
    #mddb_utils.run_sql("alter role {0} with password null".format(self.param_dict['dbuser']), dbname)
  
    for f in init_files:
      mddb_utils.load_sql("{0}/{1}".format(dbdir,f), dbname)
  
    #conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
    #                             self.param_dict['dbuser'],
    #                             self.param_dict['dbhost'],
    #                             self.param_dict['dbpass'])
  
    #mddb_utils.load_machine_info(conn, ['mddb','mddb2','qp2','qp3','qp4','qp5','qp6'])
    #conn.commit()
    #conn.close()
  
  # load stored procedure files in 'sp_files'
  # anything inside sp_files shouldn't delete existing data
  # e.g., definitions of stored procedures.
  # However, if you want to reset job control information
  # everytime you reload stored procedures then you may
  # include 'job_control.sql' here.
  def load_stored_procedures(self, sp_files):
    dbname = self.param_dict['dbname']
    dbdir = os.path.join(os.getenv('HOME'),
                         self.param_dict['user'],
                         self.param_dict['mwd'],
                         self.param_dict['dbdir'])
    print "Setting up database ", dbname
    for f in sp_files:
      mddb_utils.load_sql("{0}/{1}".format(dbdir,f), dbname)
  
  # iterates through the 'Workers' table and terminates all
  # existing worker screen sessions
  def quit_existing_workers(self):
    conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                 self.param_dict['dbuser'],
                                 self.param_dict['dbhost'],
                                 self.param_dict['dbpass'])
    cur = conn.cursor()
    cur.execute("select worker_name from Workers")
    for rec in cur:
      gw_name = rec[0]
      mddb_utils.run_cmd("screen -X -S {0} quit".format(gw_name))
  
    cur.execute("truncate Workers cascade")
    cur.execute("select setval((select pg_get_serial_sequence('workers', 'worker_id')), 1, false)")
    conn.commit()
    cur.close()
    conn.close()
  

  # start a worker screen session (each worker screen session runs a worker process) 
  # according to the parameters specified in self.param_dict
  def start_gmw(self, conn):
    self.param_dict['gmw'] = os.path.join(os.getenv('HOME'), 
                             self.param_dict['user'], 
                             self.param_dict['mwd'], 
                             self.param_dict['gmw'])

    print os.environ
    sys.stdout.flush()

    args =       """--worker_name {worker_name} """+\
                 """--worker_id {worker_id} """+\
                 """--dbpass {dbpass} """+\
                 """--dbname {dbname} --dbhost {dbhost} --dbuser {dbuser} """+\
                 """--gen_host {gen_host} """ +\
                 """--prefix {prefix}/{dbname} """+\
                 """--mwd {mwd} """ +\
                 """--user {user} """ +\
                 """--avg_seq_jobs {avg_seq_jobs} """ +\
                 """--res_config_name {res_config_name} """ 
  
    # start a screen here
    start_local_screen =  """screen -S {worker_name} -d -m """+\
                 """env PYTHONPATH={pythonpaths} GMX_MAXBACKUP=-1 GPHOME=/opt/greenplum-db/. """ +\
                 """PATH={PATH} LD_LIBRARY_PATH={LD_LIBRARY_PATH} SSH_AUTH_SOCK={SSH_AUTH_SOCK} """.format(**os.environ) +\
                 """/opt/greenplum-db/ext/python/bin/python {gmw} """ + args
  
    # ssh into the worker host and start a screen
    start_remote_screen = """screen -S {worker_name} -d -m ssh {worker_host} """+\
                 """env PYTHONPATH={pythonpaths} GMX_MAXBACKUP=-1 """ +\
                 """"/usr/bin/python {gmw} """ + args  + """" """


    template = start_local_screen if self.param_dict['worker_host'] == 'localhost' else start_remote_screen
    mddb_utils.run_cmd(template.format(**self.param_dict))
  
  # make sure that the corresponding remote resource has up-to-date scripts
  def prepare_worker(self, h):
    remote_resource_names = [k for k,v in resources.resource_dict.items() if v != resources.LocalResource]
    if h in remote_resource_names:
      res_class = resources.resource_dict[h]
      res_paths   = res_class.get_paths()
      script_dirs = res_class.script_dirs or ['mddb/scheduler','mddb/templates/']
      resources.Resource.sync_scripts(res_class.gateway_host, 
                                      res_paths['resource_prefix'], 
                                      self.param_dict['user'],
                                      script_dirs)

  # iterates through a list 'l' of workers and starts worker processes
  # there are three type of resources: localhost, localcluster-remotehost,
  # and remote-resource. For the first two, the worker screen session
  # runs on the worker host 'worker_host' and the data is also generated
  # on the woker host. For remote-resource, the worker screen session is
  # running on localhost but the data is generated on the remote resource.
  def setup_workers(self, l):
    conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                 self.param_dict['dbuser'],
                                 self.param_dict['dbhost'],
                                 self.param_dict['dbpass'])
  
    for d in l:
      h = d['hostname']
      self.param_dict['res_config_name'] = d.get('res_config_name') or 'default'
      self.param_dict['avg_seq_jobs'] = d.get('avg_seq_jobs') or 1
      self.prepare_worker(h)
      prefix = os.path.join(os.getenv('HOME'),
                            self.param_dict['user'],
                            self.param_dict['mwd'])

      self.param_dict['remotegen'] = ""
      
      self.param_dict['gen_host'] = h
      remote_resource_names = [k for k,v in resources.resource_dict.items() if v != resources.LocalResource]
      if h in remote_resource_names:
        self.param_dict['worker_host'] = 'localhost'
      else:
        self.param_dict['worker_host'] = self.param_dict['gen_host']

      # centos hosts reqires a different PYTHONPATH than ubuntu hosts
      # TODO: add more centos hosts all qp-hm* and qp-hd*
      if self.param_dict['worker_host'] in ['mddb', 'mddb2', 'localhost']:
        self.param_dict['pythonpaths'] = self.centos_pythonpaths.format(prefix)
      else:
        self.param_dict['pythonpaths'] = self.ubuntu_pythonpaths.format(prefix)

      # insert the worker name into the 'workers' table
      i = mddb_utils.execute_query(conn, "select workers_insert('{gen_host}');".format(**self.param_dict))
      self.param_dict['worker_id']  = i
      self.param_dict['worker_name'] = "gmw_{dbname}_{gen_host}_{worker_id}_{user}".format(**self.param_dict)
      self.start_gmw(conn)
    conn.commit()
    conn.close()
  
  def quit_everything(self):
    try:
      self.quit_existing_workers()
      # quit the control loop too if implemented
    except:
      pass

  # check whether the user name is the same as the shared account name
  # if not, sync the user's mddb directory inside the shared account's home directory
  # and launch the exact same command as the shared account user
  @staticmethod
  def check_user():
    print 'check_user: ', os.environ['USER'] 
    if os.environ['USER'] != Controller.sharedacc:
      cmd = 'ssh {0} mkdir {1}'.format(Controller.sharedacc, os.getenv("USER"))
      print cmd
      subprocess.call(cmd, shell = True)

      cmd = 'rsync -av -f"- .git/" {0}/mddb/ {1}@localhost:{2}/mddb/'.format(
                                                                   os.getenv("HOME"),
                                                                   Controller.sharedacc,
                                                                   os.getenv("USER"))
      print cmd
      subprocess.call(cmd, shell = True)

      prefix = os.path.join('/home',
                            Controller.sharedacc,
                            os.getenv('USER'),
                            Controller.mwd)


      pp = Controller.centos_pythonpaths.format(prefix)


      p   = os.getenv('PATH')
      fn  = os.path.join(os.getenv("USER"), os.path.relpath(os.getcwd(), os.getenv("HOME")), sys.argv[0])
      env = "PYTHONPATH={0} GPHOME=/opt/greenplum-db/. PGHOST=localhost4 PGOPTIONS=--client-min-messages=error ".format(pp)
      env = env + "PATH={PATH} LD_LIBRARY_PATH={LD_LIBRARY_PATH} ".format(**os.environ)
      #env = env + "SHAREDACC={0} ".format(Controller.sharedacc)
      # repeat the same command with appropriate environment variables 
      cmd = 'ssh {0} "source .profile; source .keychain/{1}-sh; env {2} python {3} {4} --user {5}"'.format(
             Controller.sharedacc,
             os.environ['HOSTNAME'], env, fn, " ".join(sys.argv[1:]), os.getenv("USER"))
      print cmd
      subprocess.call(cmd, shell = True)
      return False
    return True
