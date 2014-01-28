import subprocess
try:
  import psycopg2
except:
  pass

import time
import os
import math
import datetime


try:

  import parser as mddb_parser
  import generate_peptides as genpep
except:
  pass

import cStringIO
import itertools
import ast

import json

residue_dict = {
   'a':'ALA', 'r':'ARG', 'n':'ASN', 'd':'ASP',
   'c':'CYS', 'e':'GLU', 'q':'GLN', 'g':'GLY',
   'h':'HSD', #'x':'HSD', 'z': 'HSE',
   'i':'ILE', 'l':'LEU', 'k':'LYS',
   'm':'MET', 'f':'PHE', 'p':'PRO', 's':'SER',
   't':'THR', 'w':'TRP', 'y':'TYR', 'v':'VAL' }


def get_residues(protein_seq):
  residues = []
  prefix = None
  for c in protein_seq.lower():
    if c in ['*', '+']:
      prefix = c
    else:
      if prefix != None:
        residues = residues + [residue_dict[prefix+c]]
      else:
        residues = residues + [residue_dict[c]]                                                  
      prefix = None
  print protein_seq, " -> ", residues
  return residues

def convert(v):
  if v == None:
    return "NULL"
  elif isinstance(v, bool):
    return "TRUE" if v else "FALSE"
  elif isinstance(v, str):
    return "'{0}'".format(v)
  elif isinstance(v, list):
    return "ARRAY[" + ",".join(map(convert, v)) + "]"
  else:
    return str(v)



def mdq_push(conn
             ,protein_id=None
             ,protein_seq=None
             ,num_expjobs=1
             ,num_iterations=0
             ,num_steps_sim=10000
             ,user_id=1
             ,policy='fixed'
             ,simulator='charmm'
             ,config_id=1
             ,trj_fn=None
             ,psf_fn=None
             ,well_fn=None
             ,timestep_pico=None
             ,trj_save_freq=None
     ):

  terms_st = ','.join(map(lambda x: "'{0}'".format(x), protein_seq))
  d = {'num_expjobs':       num_expjobs
       ,'num_iterations':   num_iterations
       ,'expanded_terms':   None #"'{}'::text[]"
       ,'rest_of_terms':    list(protein_seq)
       ,'num_terms_left':   len(protein_seq)
       ,'num_steps_sim':    num_steps_sim
       ,'timestep_pico':    timestep_pico
       ,'trj_save_freq':    trj_save_freq
       ,'user_id':          user_id
       ,'protein_id':       protein_id
       ,'policy':           policy
       ,'simulator':        simulator
       ,'config_id':        config_id
       ,'trj_fn':           trj_fn
       ,'psf_fn':           psf_fn
       ,'well_fn':          well_fn
      }  


  if d['num_iterations'] == 0:
    d['num_steps_sim'] = None
    d['policy'] = None
    d['simulator'] = None

  cols,vals = zip(*d.items())
  cols = ','.join(map("entry.{0}".format, cols))
  vals = ','.join(map(convert, vals))

  

  st = """insert into MDQueue ({0}) values ({1})""".format(cols,vals)
  execute_query(conn, st)

def load_sql(sql_file, dbname):
  some_options = "-X -q  -1 -v ON_ERROR_STOP=0 --pset pager=off"

  run_cmd("/opt/greenplum-db/bin/psql {0} -d {1} -f {2}".format(some_options, dbname, sql_file))

def run_sql(sql_cmd, dbname):
  some_options = "-X -q  -1 -v ON_ERROR_STOP=0 --pset pager=off"
  run_cmd("""/opt/greenplum-db/bin/psql {0} -d {1} -c "{2};" """.format(some_options, 
          dbname, sql_cmd))

def run_cmd(cmd):
  print cmd
  ret= subprocess.Popen(cmd, shell=True,
                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
  print ret
  return ret




def run_cmd_with_retries(cmd, n):
  c = 0
  while True:
    c = c+1
    try:
      print c, ': ', cmd
      ret= subprocess.Popen(cmd, shell=True,
                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
      break
    except Exception, e:
      time.sleep(0.2)
      print str(e)

    if c > n:
      break;


def run_cmd_background(cmd):
  print cmd
  subprocess.Popen(cmd, shell=True,
                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


def run_cmd_with_file(cmd, fn, log_fn):
  #print log_fn
  with open(log_fn, "w") as debug_log:
    print cmd, " < ", fn
    ret = ''
    with open(fn, 'r') as myinput:
      ret= subprocess.Popen(cmd, shell=True, stdin=myinput,
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    debug_log.write(ret + "\n")
    #print ret
  return ret

def run_cmd_with_env(cmd, log_fn, my_env):
  print cmd

  with open(log_fn, "w") as debug_log:
    ret= subprocess.Popen(cmd, shell=True, env = my_env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()

    debug_log.write(ret + "\n")
  return ret

def get_dbconn(dbname, dbuser, dbhost, gp_pass):
  conn = None
  gp_pass = open(gp_pass, 'r').read()
  count = 1
  while conn == None:
    try:
      conn = psycopg2.connect(database=dbname, user=dbuser, host=dbhost,
                           password=gp_pass.decode('base64').decode('rot13'))
    except Exception, e:
      message = str(e)
      if 'not exist' in message:
        return None
      count = count + 1
      print "Trial {0}: {1}".format(count, e)
      time.sleep(1)
      conn = None

  #print "Connected to ", dbname
  return conn

def execute_query(conn, st):
  cur = conn.cursor()
  print st
  cur.execute(st)
  retval = None
  try:
    retval = cur.fetchone()[0]
  except:
    pass

  cur.close()
  conn.commit()
  return retval

def pbe_split(st):
  st = st.strip()
  p, e = os.path.splitext(st)
  p, b = os.path.split(path)
  return p, b , e

def modify_namd_inp(in_fn, trj_id):
  path, ext = os.path.splitext(in_fn)
  out_fn = path + "." + str(trj_id) + ext
  lines = None
  with open(in_fn, 'r') as in_file:
    lines = in_file.readlines()
  st = 'coordinates'
  with open(out_fn, 'w') as out_file:
    for l in lines:
      if l.startswith(st):
        path, ext = os.path.splitext(l.split()[1])
        new_pdb_fn = path + "." + str(trj_id) + ext 
        out_file.write(st + "        " + new_pdb_fn+ "\n")
      else:
        out_file.write(l)
  return out_fn, new_pdb_fn


def run_charmm(charmm_in_fn):
  charmm_bin = '/damsl/mddb/software/MD/Charmm/c36b1/exec/gnu/charmm'
  st =run_cmd_with_file(charmm_bin, charmm_in_fn, charmm_in_fn+".log")
  return st

def load_protein(conn, c_id, p_id, in_path_name, prm_format, csv_prefix, user):

  if csv_prefix == None:
    psf_dir,psf_fn = os.path.split(in_path_name)
    run_cmd("rm {0}/ch_*.csv".format(psf_dir))

    cmd = 'python /home/webapp/{0}/mddb/simulation/parser.py --mode protein'.format(user)+\
                 ' --trjtop {0} --inputdir {1} --outputdir {2} --sim {3} ch 1'.format(psf_fn, psf_dir, psf_dir, prm_format)

    st = run_cmd(cmd)
    param_prefix = "{0}/ch_{1}".format(psf_dir, prm_format)
    execute_query(conn, "update Proteins set protein_fn = '{0}' where protein_id = {1}".format(in_path_name, p_id))
  else:
    param_prefix = csv_prefix

  execute_query(conn, "select load_protein({0},{1},'{2}')".format(
                c_id, p_id, param_prefix));

def queue_gmjobs(conn, task_id, job_lists):
  cur = conn.cursor()
  for l in job_lists:
    sql_st = 'select queue_gearman_job(\'{0}\', \'{1}\')'.format(task_id, json.dumps(l))
    print sql_st
    cur.execute(sql_st)

  cur.close()
  conn.commit()


def gen_psf(param_dict):
  job_dict = genpep.generate_scripts({
    'dbpass':         param_dict['dbpass'],
    'dbhost':         param_dict['dbhost'],
    'dbuser':         param_dict['dbuser'],
    'dbname':         param_dict['dbname'],
    'prefix':         param_dict['prefix'] + "/" + param_dict['dbname'],
    'protein_seq':    param_dict.get('protein_seq') or "",
    'trj_id':         0,
    'protein_id':     param_dict['protein_id'],
    'sim_ff':         'all' if param_dict['config_id'] == 1 else 'ua',
  })

  print "creating parameter files", job_dict['psf_file_name']
  run_charmm(job_dict['minimization_in_fn'])
  return job_dict['psf_file_name']

def add_protein(conn, param_dict):
  c_id = param_dict.get('config_id') or 1
  psf_fn = param_dict.get('psf_fn') 
  protein_seq = param_dict.get('protein_seq')
  if protein_seq == None:
    p_id = execute_query(conn,
              "select insert_new_protein(Null, '{0}', Null, Null)".format(param_dict['protein_name']))
  else:
    p_id = execute_query(conn,
              "select insert_new_protein('{0}', Null, Null, Null)".format(param_dict['protein_seq']))

  param_dict['protein_id'] = p_id

  csv_prefix = param_dict.get('csv_prefix')

  if psf_fn == None and csv_prefix == None:
    psf_fn = gen_psf(param_dict)


  prm_format = "charmm"
  load_protein(conn, c_id, p_id, psf_fn, prm_format, csv_prefix, param_dict['user'])
  conn.commit()
  return p_id

def load_trajectories(conn, ss_id, param_dict):
  e_id   = execute_query(conn,
               """select expjobs_insert({0}, Null, Null, Null, Null, Null)""".format(ss_id))

  trj_list_fn = os.path.join(param_dict['trj_dir'], param_dict['trj_list_fn'])
  pdb_fn = os.path.join(param_dict['pdb_dir'], param_dict['pdb_fn'])

  with open(trj_list_fn, 'r') as ifp:
    fn_list = ifp.readlines()

  fn_list = map(lambda t: t.split(',')[1].rstrip(), fn_list)

  

  for dcd_fn in fn_list:
    trj_id = execute_query(conn,
                 """select trajectories_insert({0}, Null, Null, Null, Null)""".format(e_id))
    trj_fn = os.path.join(param_dict['trj_dir'], dcd_fn)



    step = 1
    l = mddb_parser.get_trjectory_length(pdb_fn, trj_fn)
    chunk_size = 10
    num_chunks = int(math.ceil(float(l)/chunk_size))
    cur = conn.cursor()
    for i in xrange(0, num_chunks):
      start = i * chunk_size
      stop  = min(start + chunk_size, l)

      print "trj_id: ", trj_id, "      range: ", start, "->", stop \
            ,"    total: ", l, "         time: ", datetime.datetime.now()

      st_io  = mddb_parser.parse_trajectory_range_to_buffer(trj_id, pdb_fn, trj_fn, start, stop, step)

      cur = conn.cursor()
      cur.copy_from(st_io, 'AtomPositions', columns=('trj_id', 't', 'atom_id', 'x', 'y', 'z'))
      cur.close()
      conn.commit()

  return

def load_wells(conn, ss_id, well_fn):
  with open(well_fn, 'r') as ifp:
    content = ifp.read()

  print content
  well_dict = ast.literal_eval(content)
  print well_dict

  cur = conn.cursor()
  for k,v in well_dict.items():
    if k == 'phipsi':
      st_io = cStringIO.StringIO()
      for w in v:
        st_io.write("{subspace_id},{well_id},{ref_id},{feature_id},{centroid},{span}\n".format(subspace_id=ss_id,**w))
      print st_io.getvalue()
      st_io.seek(0)
      cur.copy_from(st_io, 'MSM_EnergyWellsPhiPsi', sep=',')

  cur.close()
  conn.commit()
  return


def subspace_init(param_dict):
  conn = get_dbconn(param_dict['dbname'],
                    param_dict['dbuser'],
                    param_dict['dbhost'],
                    param_dict['dbpass'])

  if param_dict.get('protein_id') in [None, 0]:
    param_dict['protein_id'] = add_protein(conn, param_dict)

  ss_id = execute_query(conn, 
                        """select subspaces_insert({config_id}, {protein_id}, '{gmdhost}', 
                                                   Null, '{user_id}', '{policy}');""".format(**param_dict))
  print "subspace_id: ", ss_id
  param_dict['subspace_id'] = ss_id


  if param_dict.get('trj_list_fn') != None:
    load_trajectories(conn, ss_id, param_dict)

  if param_dict.get('well_fn') != None:
    load_wells(conn, ss_id, param_dict['well_fn'])



  if param_dict.get('num_expjobs') != None:
    for _ in itertools.repeat(None, param_dict['num_expjobs']):
      e_id   = execute_query(conn,
                """select expjobs_insert({subspace_id},{num_steps_sim},{timestep_pico},
                                         {trj_save_freq},{num_iterations},'{simulator}')""".format(**param_dict))
    execute_query(conn, "select subspaces_activate({0})".format(ss_id))

  conn.close()
  return ss_id




def latest_file(dir):
  files = [os.path.join(dir, fname) for fname in os.listdir(dir)]
  latest = max(files, key=os.path.getmtime)
  return latest

def read_dict(st):
  dict = {}
  for line in st.split("\n"):
    terms = line.split(":")
    if len(terms) != 2:
      continue
    dict[terms[0]] = terms[1].strip()

  return dict

def load_machine_info(conn, machine_list):
  cur = conn.cursor()
  sql = "truncate MachineInfo"
  cur.execute(sql)
  for m in machine_list:
    st = run_cmd("ssh {0} lscpu".format(m))
    st = st + run_cmd("ssh {0} cat /proc/meminfo | grep MemTotal".format(m))
    dict = read_dict(st)
    for k,v in dict.iteritems():
      sql = "insert into MachineInfo values (E'{0}', E'{1}',E'{2}')".format(m, k,v)
      cur.execute(sql)
  conn.commit()

def check_dcd(trj_id, psf_fn, dcd_fn, expected_dcd_length):
    progress = 0
    count = 0
    while progress < expected_dcd_length:
      time.sleep(10)
      try:
        u = MDAnalysis.Universe(psf_fn, dcd_fn)
        progress = len(u.trajectory)
        print progress, " frames generated"
        if count % 20 == 0:
          conn = mddb_utils.get_dbconn(dbname, dbuser, dbhost, dbpass)
          update_benchmark_entry(conn, trj_id, progress)
          conn.close()
        count = count+1
      except:
        pass


def gen_conversion_script(in_pdb_fn, out_pdb_fn, script_fn):

  st = """library(bio3d)\n""" +\
       """pdb <- read.pdb('{0}')\n""".format(in_pdb_fn) +\
       """new <- convert.pdb(pdb, type="amber")\n""" +\
       """write.pdb(new, file='{0}')\n""".format(out_pdb_fn)

  open(script_fn, 'w').write(st)

def reformat_pdb1(pdb_fn):
  #mddb_utils.run_cmd("/damsl/mddb/software/tools/MMTSB/mmtsb/perl/convpdb.pl -out amber " +\
  #        "{0} > {0}2".format(pdb_fn))
  mddb_utils.run_cmd("cp {0} {0}.bak".format(pdb_fn))
  script_fn =  pdb_fn + '.R'
  gen_conversion_script(pdb_fn, pdb_fn + '2', script_fn)
  mddb_utils.run_cmd("Rscript " + script_fn)
  mddb_utils.run_cmd("mv {0}2 {0}".format(pdb_fn))


