import os, os.path, subprocess, tempfile
import datetime, time
from optparse import OptionParser
import psycopg2
import sys
import parser as mddb_parser
import json
import base64
import xml.etree.ElementTree as ET
import cStringIO
import socket
import MDAnalysis
import mddb_utils
import traceback
from threading import Thread
import generate_peptides as genpep
import re
import math
import random
import lddb_parsers
import shutil
import platform
import qsub
import gearman

omega2_bin = '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/omega2'
dock_bin   = '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/dock6'
prepare_amber_pl = '/damsl/projects/molecules/software/DD/Dock/Dock/dock6/bin/prepare_amber.pl'
in_dir     = '/damsl/projects/molecules/data/GCPII/oe/ligands/sdf_lib4'
out_dir    = '/damsl/projects/molecules/data/GCPII/oe/ligands/sdf_lib4'

fred_bin   = '/damsl/projects/molecules/software/DD/OpenEye/centos/openeye/bin/fred'
if 'Ubuntu' in platform.platform():
  fred_bin = '/damsl/projects/molecules/software/DD/OpenEye/openeye/bin/fred'

dbname = None
dbhost = None
dbpass = None
dbuser = None
gmdhost = None
gmdport = None

amber_in = """
  ligand_atom_file                                             {mol_fn}
  limit_max_ligands                                            no
  skip_molecule                                                no
  read_mol_solvation                                           no
  calculate_rmsd                                               no
  use_database_filter                                          no
  orient_ligand                                                no
  use_internal_energy                                          no
  flexible_ligand                                              no
  bump_filter                                                  no
  score_molecules                                              yes
  contact_score_primary                                        no
  contact_score_secondary                                      no
  grid_score_primary                                           no
  grid_score_secondary                                         no
  dock3.5_score_primary                                        no
  dock3.5_score_secondary                                      no
  continuous_score_primary                                     no
  continuous_score_secondary                                   no
  descriptor_score_primary                                     no
  descriptor_score_secondary                                   no
  gbsa_zou_score_primary                                       no
  gbsa_zou_score_secondary                                     no
  gbsa_hawkins_score_primary                                   no
  gbsa_hawkins_score_secondary                                 no
  amber_score_primary                                          yes
  amber_score_secondary                                        no
  amber_score_receptor_file_prefix                             3D7F_h
  amber_score_movable_region                                   ligand
  amber_score_minimization_rmsgrad                             0.01
  amber_score_before_md_minimization_cycles                    {min_steps}
  amber_score_md_steps                                         {md_steps}
  amber_score_after_md_minimization_cycles                     {min_steps}
  amber_score_gb_model                                         5
  amber_score_nonbonded_cutoff                                 18.0
  amber_score_temperature                                      300.0
  amber_score_abort_on_unprepped_ligand                        yes
  ligand_outfile_prefix                                        output
  write_orientations                                           yes
  num_scored_conformers                                        1
  rank_ligands                                                 no
"""

dock6_in = """
  ligand_atom_file                                             {mol_fn}
  limit_max_ligands                                            no
  skip_molecule                                                no
  read_mol_solvation                                           no
  calculate_rmsd                                               no
  use_database_filter                                          no
  orient_ligand                                                yes
  automated_matching                                           yes
  receptor_site_file                                           selected_spheres.sph
  max_orientations                                             500
  critical_points                                              no
  chemical_matching                                            no
  use_ligand_spheres                                           no
  use_internal_energy                                          yes
  internal_energy_rep_exp                                      12
  flexible_ligand                                              no
  user_specified_anchor                                        no
  limit_max_anchors                                            no
  min_anchor_size                                              40
  pruning_use_clustering                                       yes
  pruning_max_orients                                          100
  pruning_clustering_cutoff                                    100
  pruning_conformer_score_cutoff                               25.0
  use_clash_overlap                                            no
  write_growth_tree                                            no
  write_orientations                                           yes
  bump_filter                                                  no
  score_molecules                                              yes
  contact_score_primary                                        no
  contact_score_secondary                                      no
  grid_score_primary                                           yes
  grid_score_secondary                                         no
  grid_score_rep_rad_scale                                     1
  grid_score_vdw_scale                                         1
  grid_score_es_scale                                          1
  grid_score_grid_prefix                                       grid
  dock3.5_score_secondary                                      no
  continuous_score_secondary                                   no
  descriptor_score_secondary                                   no
  gbsa_zou_score_secondary                                     no
  gbsa_hawkins_score_secondary                                 no
  amber_score_secondary                                        no
  minimize_ligand                                              yes
  minimize_anchor                                              yes
  minimize_flexible_growth                                     yes
  use_advanced_simplex_parameters                              no
  simplex_max_cycles                                           1
  simplex_max_iterations                                       1000
  simplex_tors_premin_iterations                               0
  simplex_score_converge                                       0.1
  simplex_cycle_converge                                       1.0
  simplex_trans_step                                           1.0
  simplex_rot_step                                             0.1
  simplex_tors_step                                            10.0
  simplex_anchor_max_iterations                                500
  simplex_grow_max_iterations                                  500
  simplex_grow_tors_premin_iterations                          0
  simplex_random_seed                                          0
  simplex_restraint_min                                        no
  atom_model                                                   all
  vdw_defn_file                                                /damsl/projects/molecules/software/DD/Dock/Dock/dock6/parameters/vdw_AMBER_parm99.defn
  flex_defn_file                                               /damsl/projects/molecules/software/DD/Dock/Dock/dock6/parameters/flex.defn
  flex_drive_file                                              /damsl/projects/molecules/software/DD/Dock/Dock/dock6/parameters/flex_drive.tbl
  ligand_outfile_prefix                                        anchor_and_grow
  write_orientations                                           no
  num_scored_conformers                                        1
  rank_ligands                                                 no"""

def compose_dock6_in(in_dict, fn):
  out_content = None
  if in_dict['runamber']:
    out_content = amber_in.format(**in_dict)
  else:
    out_content = dock6_in.format(**in_dict)

  with open(fn, "w") as ofs:
    ofs.write(out_content)

def run_dock6_with_amber(mol_fn, pdb_fn, dock6_inp, dock6_out, sub_dir):
  # Example: prepare_amber.pl 70_charged.mol2 3D7F_h.pdb
  #          dock6 -i dock_70_2nd_b.in -v
  bname,ext = os.path.splitext(mol_fn)
  mol_amber_score_fn = "{0}.amber_score{1}".format(bname, ext)
  cmd_st = "cp {0} {1}".format(mol_fn, mol_amber_score_fn)
  print cmd_st
  ret = subprocess.Popen(cmd_st, cwd = in_dir + "/" + sub_dir,
                         shell=True, stdout=subprocess.PIPE,
                         stdin=subprocess.PIPE).stdout.read()

  cmd_st = 'echo "@<TRIPOS>AMBER_SCORE_ID\n{0}" >> {1}'.format(
             bname, mol_amber_score_fn)
  print cmd_st
  ret = subprocess.Popen(cmd_st, cwd = in_dir + "/" + sub_dir,
                         shell=True, stdout=subprocess.PIPE,
                         stdin=subprocess.PIPE).stdout.read()
  inp_dict = {"runamber": True,
              "mol_fn": mol_amber_score_fn,
              "recp_name": os.path.splitext(pdb_fn)[0],
              "min_steps": 500,
              "md_steps": 30000,
              }

  compose_dock6_in(inp_dict, dock6_inp)



  print ret
  cmd_st = "nice {0} -i {1} -v > {2}".format(dock_bin, dock6_inp, dock6_out)
  print cmd_st
  ret = subprocess.Popen(cmd_st, cwd = in_dir + "/" + sub_dir,
                         shell=True, stdout=subprocess.PIPE,
                         stdin=subprocess.PIPE).stdout.read()
  print ret

def run_dock6(mol_fn, dock6_inp, dock6_out, sub_dir):
  # dock6 -i dock.in > dock.out &
  inp_dict = {"runamber": False,
              "mol_fn": mol_fn,
              "recp_sph": "selected_spheres.sph",
              "receptor_id": receptor_id}

  compose_dock6_in(inp_dict, dock6_inp)

  cmd_st = "nice {0} -i {1} -v > {2}".format(dock_bin, dock6_inp, dock6_out)
  print cmd_st
  ret = subprocess.Popen(cmd_st, cwd = in_dir + "/" + sub_dir,
                         shell=True, stdout=subprocess.PIPE,
                         stdin=subprocess.PIPE).stdout.read()
  print ret

def run_fred(recp_fn, mol_fn, fred_prefix, sub_dir):
  cmd_st = ("nice {0} -receptor {1} -dbase  {2} " +
            "-save_component_scores -dock_resolution High -prefix {3}"
           ).format(fred_bin, recp_fn, mol_fn, fred_prefix)
  print cmd_st
  subprocess.Popen(cmd_st, cwd = in_dir,
                   shell=True, stdout=subprocess.PIPE,
                   stdin=subprocess.PIPE).stdout.read()




def insert_fred_scores(r_id, c_id, p_id, chargeset, tautomers, scores):
  if scores != None and len(scores) != 0:
    vals = '({0})'.format(', '.join(map(str, [r_id,c_id, p_id, chargeset, tautomers] + scores))) 
    stm = "insert into BindingScoresFred values\n{0}".format(vals)
  else:
    stm = "insert into BindingScoresFred (receptor_id, compound_id, pose_id, chargeset_id, tautomers) " +\
          "values ({0},{1},{2},{3},{4})".format(r_id,c_id, p_id, chargeset, tautomers)

  conn = mddb_utils.get_dbconn(dbname, dbuser, dbhost, dbpass)
  mddb_utils.execute_query(conn, stm)



def insert_dock6_scores(r_id, c_id, p_id, chargeset, tautomers, list_of_score_tups):

  if list_of_score_tups != None and len(list_of_score_tups) != 0:
  
    vals = ',\n'.join(
      map(lambda x: '({0})'.format(', '.join(map(str, [r_id,c_id, p_id, chargeset, tautomers] + list(x)))), list_of_score_tups)
    )
  
    stm = "insert into BindingScoresDock6 values\n{0}".format(vals)
  else:
    stm = "insert into BindingScoresDock6 (receptor_id, compound_id, pose_id, chargeset_id, tautomers) " +\
          "values ({0},{1},{2},{3},{4})".format(r_id,c_id, p_id, chargeset, tautomers)


  conn = mddb_utils.get_dbconn(dbname, dbuser, dbhost, dbpass)
  mddb_utils.execute_query(conn, stm)

def insert_dock6_amber_scores(r_id, c_id, p_id, chargeset, tautomers, list_of_score_tups):
  if list_of_score_tups != None and len(list_of_score_tups) != 0:
    vals = ',\n'.join(
      map(lambda x: '({0})'.format(', '.join(map(str, [r_id,c_id, p_id, chargeset, tautomers] + list(x)))), list_of_score_tups)
    )
  
    stm = "insert into BindingScoresDock6Amber values\n{0}".format(vals)
  else:
    stm = "insert into BindingScoresDock6Amber (receptor_id, compound_id, pose_id, chargeset_id, tautomers) " +\
          "values ({0},{1},{2},{3},{4})".format(r_id,c_id, p_id, chargeset, tautomers)

  conn = mddb_utils.get_dbconn(dbname, dbuser, dbhost, dbpass)
  mddb_utils.execute_query(conn, stm)


def prepare_filenames(in_dir, sub_dir, c_id, p_id):
  fn_dict = {
    'ambr_mol_fn':       "{0}/{1}/{2}_amberff94_{3}.mol2".format(in_dir, sub_dir, c_id, p_id),
    'gast_mol_fn':       "{0}/{1}/{2}_gasteiger_{3}.mol2".format(in_dir, sub_dir, c_id, p_id),
    'mmff_mol_fn':       "{0}/{1}/{2}_mmff94_{3}.mol2".format(in_dir,  sub_dir,c_id, p_id),
    'taut_ambr_mol_fn':  "{0}/{1}/{2}_tautomers_amberff94_{3}.mol2".format(in_dir, sub_dir, c_id, p_id),
    'taut_gast_mol_fn':  "{0}/{1}/{2}_tautomers_gasteiger_{3}.mol2".format(in_dir, sub_dir, c_id, p_id),
    'taut_mmff_mol_fn':  "{0}/{1}/{2}_tautomers_mmff94_{3}.mol2".format(in_dir, sub_dir, c_id, p_id),
    'dock6_inp':         "{0}/{1}/{2}_dock6_{3}.in".format(in_dir, sub_dir, c_id, p_id),
    'dock6_out':         "{0}/{1}/{2}_dock6_{3}.out".format(in_dir, sub_dir, c_id, p_id),
    'dock6_amber_inp':   "{0}/{1}/{2}_dock6_amber_{3}.in".format(in_dir, sub_dir, c_id, p_id),
    'dock6_amber_out':   "{0}/{1}/{2}_dock6_amber_{3}.out".format(in_dir, sub_dir, c_id, p_id)
  }
  return fn_dict;

def prepare_filenames(job_dict):
  prefix = "{in_dir}/{c_id}/r{r_id}_c{c_id}_p{p_id}_cs{cs_id}_t{taut}_".format(**job_dict)
  fn_dict = {
    'dock6_inp':         "{0}_dock6.in".format(prefix),
    'dock6_out':         "{0}_dock6.out".format(prefix),
    'dock6_amber_inp':   "{0}_dock6_amber.in".format(prefix),
    'dock6_amber_out':   "{0}_dock6_amber.out".format(prefix),
    'fred_prefix':       "{0}_fred".format(prefix),
    'fred_out':          "{0}_fred_score.txt".format(prefix),
    'qsub':              "{0}_qsub.txt".format(prefix),
    'job_json':          "{0}_job_json.txt".format(prefix)
  }
  return fn_dict



def result_exists(conn, tn, r_id, c_id, p_id, chargeset, tautomers):
  stm = "select 1 from {0} where receptor_id = {1} and compound_id = {2} and pose_id = {3} and chargeset_id = {4} and tautomers = {5}".format(tn, r_id, c_id, p_id, chargeset, tautomers)
  return mddb_utils.execute_query(conn, stm) != None

def num_finished_remote_jobs():
  n = mddb_utils.run_cmd("""ssh hhpc 'sqlite3 /home/yahmad/lddb_jobs/db/docking.db  "select count(*) from FinishedJobsView"' """)
  return int(n)

def fetch_results():
    st = mddb_utils.run_cmd("ssh hhpc python2.6 nutanong/mddb/scheduler/lddb_sqlite_utils.py --mode getresults")
    sql_stms = st.split("\n")

    conn = mddb_utils.get_dbconn(dbname, dbuser, dbhost, dbpass)
    cur = conn.cursor()
    for stm in sql_stms:
      print stm
      try:
        cur.execute(stm)
      except:
        pass

      conn.commit()
    cur.close()
    conn.close()



def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def queue_more_jobs(job_dict):
  print "queuing more jobs"
  offset = job_dict['offset']
  limit  = job_dict['limit']
  docker = job_dict['docker']

  stm = 'select receptor_id, model_id, fred_fn, dock6_fn, compound_id, tautomers, smi_fn, sdf_fn ' +\
        'from UnfinishedJobsView J ' +\
        'order by job_id ' +\
        'offset {0} limit {1}'.format(offset, limit)

  conn = mddb_utils.get_dbconn(dbname, dbuser, dbhost, dbpass)
  cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
  cur.execute(stm)
  rows = cur.fetchall()

  job_size = 1
  num_entries = len(rows)
  groups = list(chunks(rows, job_size))
  for rows in groups:

    json_st = "[{0}]".format(",".join(rows))

    stm = "select queue_gearman_job(E'{0}', E'{1}', E'{2}')".format(
                           json_st,
                           gmdhost, 
                           gmdport)

    print stm
    cur.execute(stm)


  cur.close()
  conn.commit()
  conn.close()  

def get_amber_file_names(compound_dir, receptor_name, mol_fn):
  mol_fn_rel = os.path.basename(mol_fn)
  bname,m2ext = os.path.splitext(mol_fn_rel)
  suffs1 = ['mol2', 'amber_score.mol2', 'amber.pdb', 'frcmod', 
           'gaff.mol2', 'inpcrd', 'log', 'mopac.out', 'prmtop']
  suffs2 = ['amber.pdb', 'inpcrd', 'log', 'prmtop']
  l = map(lambda x: compound_dir + "/" + bname + "." + x, suffs1) +\
      map(lambda x: compound_dir + "/" + receptor_name +"." + bname + "." + x, suffs2)
  return l


def get_chargeset_id(mol_fn):
  chargeset = None
  if 'gasteiger' in mol_fn:
    chargeset = 1
  elif 'amberff94' in mol_fn:
    chargeset = 2
  elif 'mmff94' in mol_fn:
    chargeset = 3
  return chargeset

def run_docking_combo(job_dict):
  #mol_fn     = in_dir + "/" + job_dict['path']
  mol_fn     = job_dict['path']

  taut  = 'tautomer' in mol_fn
  cs_id = get_chargeset_id(mol_fn)

  r_id       = 1

  job_dict['in_dir'] = in_dir
  job_dict['r_id'] = r_id
  job_dict['cs_id']= cs_id
  job_dict['taut'] = taut
  job_dict['mol_fn'] = mol_fn

  mol_fn_rel         = os.path.basename(mol_fn)
  job_dict['mol_fn_rel'] = mol_fn_rel

  if remotegen is None:
    run_docking_combo_local(job_dict)
  else:
    run_docking_combo_remote(job_dict)

def run_docking_combo_remote(job_dict):
  fn_dict = prepare_filenames(job_dict)
  job_json_fn = fn_dict['job_json']
  remote_dir = mddb_utils.run_cmd("ssh hhpc 'echo -n $INDIR'")
  print job_dict

  qsub.create_lddb_job_json(job_json_fn, job_dict)
  job_json_fn_remote = remote_dir + "/" + os.path.basename(job_json_fn)
  qsub_fn     = fn_dict['qsub']
  qsub.create_lddb_qsub_script(qsub_fn, job_json_fn_remote)
  qsub_fn_remote = remote_dir + "/" + os.path.basename(qsub_fn)


  mol_fn  = job_dict['mol_fn']
  sub_dir = str(job_dict['c_id'])

  fl = get_amber_file_names(in_dir + "/" + sub_dir, receptor_id, mol_fn)
  
  fl = fl + [job_json_fn, qsub_fn]
  map(lambda x: mddb_utils.run_cmd("scp {0} hhpc:{1}".format(x, remote_dir)), fl)
  mddb_utils.run_cmd("ssh hhpc qsub {0}".format(qsub_fn_remote))
  time.sleep(50)
  return ''


def run_docking_combo_local(job_dict):
  r_id               = job_dict['r_id']
  c_id               = job_dict['c_id']
  p_id               = job_dict['p_id']
  cs_id              = job_dict['cs_id']
  taut               = job_dict['taut'] 
  mol_fn             = job_dict['mol_fn']
  docker             = job_dict['docker']

  mol_fn_rel         = os.path.basename(mol_fn)


  sub_dir            = "{0}".format(c_id)
  if not os.path.exists(in_dir + "/" + sub_dir + "/grid.bmp"):
    shutil.copy(in_dir + "/grid.bmp", in_dir + "/" + sub_dir)

  if not os.path.exists(in_dir + "/" + sub_dir + "/grid.nrg"):
    shutil.copy(in_dir + "/grid.nrg", in_dir + "/" + sub_dir)

  if not os.path.exists(in_dir + "/" + sub_dir + "/selected_spheres.sph"):
    shutil.copy(in_dir + "/selected_spheres.sph", in_dir + "/" + sub_dir)

  head,tail          = os.path.splitext(mol_fn)
  amber_score_mol_fn = head + ".amber_score" + tail


  pdb_fn             = "{0}/3D7F_h.pdb".format(in_dir)
  recp_fn            = "{0}/3D7F.oeb.gz".format(in_dir)

  fn_dict = prepare_filenames(job_dict)
  dock6_inp = fn_dict['dock6_inp']
  dock6_out = fn_dict['dock6_out']
  dock6_amber_inp = fn_dict['dock6_amber_inp']
  dock6_amber_out = fn_dict['dock6_amber_out']
  fred_prefix = fn_dict['fred_prefix']
  fred_out = fn_dict['fred_out']

  conn = mddb_utils.get_dbconn(dbname, dbuser, dbhost, dbpass)
 
  if docker == 1:
    if not result_exists(conn, "bindingscoresfred", r_id, c_id, p_id, cs_id, taut):
      print "==========================================================================="
      print "Running fred"
  
      run_fred(recp_fn, mol_fn, fred_prefix, sub_dir)
      fred_scores = lddb_parsers.parse_fred(fred_out)
      insert_fred_scores(r_id, c_id, p_id, cs_id, taut, fred_scores)
    else:
      print "Skip fred for {0}".format(mol_fn)
  
  
  if docker == 2:
    if not result_exists(conn, "bindingscoresdock6", r_id, c_id, p_id, cs_id, taut):
      print "==========================================================================="
      print "Running dock6"
      run_dock6(mol_fn_rel, dock6_inp, dock6_out, sub_dir)
      dock6_scores = lddb_parsers.parse_dock6(dock6_out)
      insert_dock6_scores(r_id, c_id, p_id, cs_id, taut, dock6_scores)
    else:
      print "Skip dock6 for {0}".format(mol_fn)
  if docker == 3: 
    if not result_exists(conn, "bindingscoresdock6amber", r_id, c_id, p_id, cs_id, taut):
      print "==========================================================================="
      print "Running dock6 with amber"
      run_dock6_with_amber(mol_fn_rel, pdb_fn, dock6_amber_inp, dock6_amber_out, sub_dir)
      dock6_amber_scores = lddb_parsers.parse_dock6_amber(dock6_amber_out)
      insert_dock6_amber_scores(r_id, c_id, p_id, cs_id, taut, dock6_amber_scores)
    else:
      print "Skip dock6amber for {0}".format(mol_fn)
  return ''



class CustomGearmanWorker(gearman.GearmanWorker):
    def set_stop(self):
        self.stop = True

    def on_job_execute(self, current_job):
        print "Job started"
        return super(CustomGearmanWorker, self).on_job_execute(current_job)

    def on_job_exception(self, current_job, exc_info):
        print "Job failed, CAN stop last gasp GEARMAN_COMMAND_WORK_FAIL", exc_info
        traceback.print_exc()
        return super(CustomGearmanWorker, self).on_job_exception(current_job, exc_info)

    def on_job_complete(self, current_job, job_result):
        print "Job completed"
        return super(CustomGearmanWorker, self).send_job_complete(current_job, job_result)

    def after_poll(self, any_activity):
        if 'stop' in vars(self).keys():
          return False
        else:
          return True


if __name__ == '__main__':
  print "Starting up a gearman worker."

  usage = "Usage: %prog [options] <config file>"
  parser = OptionParser(usage=usage)

  parser.add_option("-c", "--client", type="string", dest="client_id",
                    default="binding_worker",
                    help="set the worker client id", metavar="#CLIENTID")

  parser.add_option("-t", "--task", type="string", dest="task_id",
                    default='lddb_worker',
                    help="set the worker task id", metavar="#TASKID")

  parser.add_option("-n", "--name", type="string", dest="worker_name",
                    default='lddb_workeri_00',
                    help="set the worker name", metavar="#NAME")

  parser.add_option("-p", "--dbpass", type="string", dest="dbpass",
                    default=None,
                    help="set the password", metavar="#GPPASS")

  parser.add_option("-g", "--gpu", action="store_true", dest="run_on_gpu",
                    default=False,
                    help="gpu", metavar="#GPU")

  parser.add_option("-d", "--dbname", type="string", dest="dbname",
                    default='lddb',
                    help="set the database", metavar="#LDDB_DBNAME")

  parser.add_option("--dbhost", type="string", dest="dbhost",
                    default='mddb',
                    help="set the host", metavar="#HOST")

  parser.add_option("--remotegen", type="string", dest="remotegen",
                    default=None,
                    help="set the remote data generator", metavar="#REMOTEGEN")

  parser.add_option("--gmdhost", type="string", dest="gmdhost",
                    default='mddb',
                    help="set the gearman host", metavar="#GMDHOST")

  parser.add_option("--gmdport", type="string", dest="gmdport",
                    default=60000,
                    help="set the gearman port", metavar="#GMDPORT")

  parser.add_option("-u", "--dbuser", type="string", dest="dbuser",
                    default=os.environ['USER'],
                    help="set the user name", metavar="#USER")

  (options, args) = parser.parse_args()

  receptor_id = '3D7F_h'

  dbname      = options.dbname
  dbhost      = options.dbhost
  dbuser      = options.dbuser
  dbpass      = options.dbpass

  remotegen   = options.remotegen

  gmdhost     = options.gmdhost
  gmdport     = options.gmdport

  run_on_gpu  = options.run_on_gpu
  task_id     = options.task_id
  worker_name = options.worker_name

  worker      = CustomGearmanWorker(["{0}:{1}".format(gmdhost,gmdport)])
  worker.set_client_id(options.client_id)
  worker.register_task(task_id, process_work)

  hostname    = socket.gethostname()
  conn = mddb_utils.get_dbconn(dbname, dbuser, dbhost, dbpass)

  cur = conn.cursor()
  sys.stdout.flush()
  #cur.execute(st)
  cur.close()
  conn.commit()
  conn.close()
  worker.work()
