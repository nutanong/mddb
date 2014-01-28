import traceback
import copy
import os
import mddb_utils
import mddb_controller
import sys
import psycopg2
import json



#  runs: pdbs x simulation scripts
#



class BDSLSS_Controller(mddb_controller.MDDBController):

  def run(self):
    mode_dict = {'load': self.load,
                 'load_trjs_dir': self.load_trjs_dir,
                 'self': self.setup,
                }

    mode = self.param_dict['mode']
    mode_dict[mode]()

  def load(self):
    output = '/damsl/projects/mddb_bdslss3/trajectories/'
    pdict = {
     'gmw':           'scheduler/mddb_worker.py'
     ,'jobfile':       ''
     ,'protein_name': 'PTEN'
     ,'config_id': 1
     ,'gmdhost':'localhost'
     ,'dbname':'mddb_bdslss3'
     ,'dbhost':'mddb'
     ,'dbuser':'jbw'
     ,'dbpass':'/home/jbw/pass.out'
     ,'user_id': '1'
     ,'trj_prefix':output
     ,'policy':'long_running' 
     ,'trj_dir': ''
    }
    self.param_dict.update(pdict)
    while (self.load_one()):
      print "Moving to next sequence"
 
  def load_trjs_dir(self):    
    conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                 self.param_dict['dbuser'],
                                 self.param_dict['dbhost'],
                                 self.param_dict['dbpass'])

    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

    cur.execute("select S.pdb_id from seq_pdb S left join trajectories T on S.pdb_id = T.pdb_id where T.pdb_id is null;")
    run_params = cur.fetchall()
    cur.close()
    conn.close()

    wanted_list = map(lambda d: int(d['pdb_id']), run_params)

    print wanted_list
    trj_dir = '/damsl/projects/{dbname}/trajectories'.format(**self.param_dict)
    pdb_id_dirs = os.listdir(trj_dir)
    #print pdb_id_dirs

    for pdb_id_dir in pdb_id_dirs:
      try:
        pdb_id_dir = int(pdb_id_dir)
        if pdb_id_dir in wanted_list:
          d = {'pdb_id': pdb_id_dir, 'trj_prefix': trj_dir}
          d.update(self.param_dict)
          self.load_trj(d)
        else:
          print 'pdb_id_dir:', pdb_id_dir, 'is already loaded'
      except Exception as e:
        print e
        pass

  def setup(self):
    self.quit_everything()
    l = []
    l.append({'hostname':'lonestar', 
              'avg_seq_jobs': 1, 
              'res_config_name': 'default'})

    self.setup_workers(l)

  # The 1% job creator!
  def createjobs(self):
    conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                 self.param_dict['dbuser'],
                                 self.param_dict['dbhost'],
                                 self.param_dict['dbpass'])

    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    cur.execute("select S.pdb_id from seq_pdb S left join trajectories T on S.pdb_id = T.pdb_id where T.pdb_id is null;")
    run_params = cur.fetchall()
    

    cur.execute("truncate table jobqueue")
    #cur.execute("truncate table loadingqueue")

    for d in run_params:
      d['generator'] = 'impsolv'
      d['dbname']    = self.param_dict['dbname']
      d['dbuser']    = self.param_dict['dbuser']
      d['dbhost']    = self.param_dict['dbhost']
      d['dbpass']    = self.param_dict['dbpass']
      data = 'ARRAY' + str([json.dumps(d)])
      cur.execute('select jobqueue_insert({0})'.format(data))
      #print d

    cur.close()
    conn.commit()
    conn.close()

  
    

  def config_ss(self):
    self.load_configs()
    print "Successfully Configured Subspace"

  def load_one(self):
    try:
      query = "SELECT jq_entry_id,data FROM loadingqueue ORDER BY jq_entry_id limit 1"

      conn = mddb_utils.get_dbconn(self.param_dict['dbname'],
                                 self.param_dict['dbuser'],
                                 self.param_dict['dbhost'],
                                 self.param_dict['dbpass'])
      cur = conn.cursor()
      cur.execute(query)
  
      row = cur.fetchone()
      id = row[0]
      
      data = row[1]
      d = eval(data)


      self.load_trj(d)   
      print "Loaded: id"
      cur.execute("DELETE FROM loadingqueue WHERE jq_entry_id = %s",(id,))

      cur.close()
       
      conn.commit() 
      conn.close() 
      return True
    except Exception as e:
      traceback.print_exc()
      print "failed" 
      return False
       
  def update_dict(self,d,seqid):
    seqid = str(seqid)
    output = os.path.join(self.param_dict['trj_prefix'],seqid)
    csv = os.path.join(output,"dcds.csv")
    f = open(csv,"w+")
    f.write("0," + "run.dcd")
    f.close()
      
    d['trj_list_fn'] =   'dcds.csv'
    d['trj_dir'] =  output
    d['pdb_dir'] =  output
    d['pdb_fn'] =  'step1_pdbreader.psf'
    d['psf_fn'] = os.path.join(output,'step1_pdbreader.psf')

  def load_trj(self,d):
    print d
    d2 = copy.deepcopy(self.param_dict)
    d2.update(d)
    d = d2
    print d
    pdbid = d['pdb_id']
    query = "SELECT SS.subspace_id, SS.seqid FROM seq_pdb SP, seq_subspace SS WHERE SP.pdb_id = %s AND SP.seqid = SS.seqid"
   
    conn = mddb_utils.get_dbconn(d['dbname'],
                                 d['dbuser'],
                                 d['dbhost'],
                                 d['dbpass'])
    cur = conn.cursor()
    cur.execute(query,(pdbid,))
    
    rows = cur.fetchall()
    ss_id = None
    seq_id = None
    if (len(rows) == 0):
      print "Creating New Subspace:"
      query = "SELECT seqid from seq_pdb where pdb_id = %s"
      cur = conn.cursor()
      cur.execute(query,(pdbid,))
      seqid = cur.fetchall()[0][0]
      cur.close()
      self.update_dict(d,pdbid)
      self.config_ss()
      ss_id = mddb_utils.subspace_init(d)
    else:
      ss_id = rows[0][0]
      seqid = rows[0][1]
      self.update_dict(d,pdbid)
      mddb_utils.load_trajectories(conn,ss_id,d)

    cur.close()
    cur = conn.cursor()
    query = "INSERT INTO seq_subspace (seqid,subspace_id) VALUES (%s, %s)"
    print query
    cur.execute(query, (seqid, ss_id))
    conn.commit()  
    print cur.fetchone()[0]

    cur.close()
    conn.close()

    print "Successfuly Loaded Trajectory"
  


if __name__ == '__main__':
  pdict = {'gmw':'scheduler/worker.py'}
  if BDSLSS_Controller.check_user():
    controller = BDSLSS_Controller(pdict)
a   controller.run()

