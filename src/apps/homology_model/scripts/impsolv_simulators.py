import os
import simulators
import impsolv_pipeline 
import shutil
import json

# extends abstractSimulator for next step

class ImpSolvSimulator(simulators.AbstractSimulator):

  def run(self, input_params):
    impsolv = impsolv_pipeline.ImpSolvPipeline(input_params)
    impsolv.run()
    with open(os.path.join(input_params['output_prefix'], 'status.txt'), 'w') as ofp:
      ofp.write('done')

  @staticmethod
  def get_output_fns(d):
    return {'status_fn': 'status.txt'}
 
  @staticmethod
  def get_sync_info(d):
    return []

  def load(self, conn, result_dir, d):
    if d.get('sim_run_id'):
      dest_dir = os.path.join(d['dest_dir'], 'trajectories', str(d['sim_run_id']))
    else:#KTWexcept:
      dest_dir = os.path.join(d['dest_dir'], 'trajectories', str(d['pdb_id']))

    if os.path.isdir(dest_dir):
      shutil.rmtree(dest_dir)
    
    shutil.copytree(result_dir, dest_dir)
    cur = conn.cursor()
    sql_st = "insert into loadingqueue (jq_entry_id, data) values ({0},'{1}')".format(d['jq_entry_id'],json.dumps(d))
    cur.execute(sql_st)
    
    cur.close() 
    conn.commit()       
    return

  @staticmethod
  def preprocess(d, local_paths):
    impsolv_pipeline.ImpSolvPipeline.preprocess(d, local_paths)

