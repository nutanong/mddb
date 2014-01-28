from abc import ABCMeta, abstractmethod, abstractproperty
import os
import cStringIO
import pipeline
import mddb_utils
import subprocess
import time
import sys
import shutil
try:
  import ModellerTools as MT
except:
  pass

class ImpSolvPipeline(pipeline.Pipeline):
  __metaclass__ = ABCMeta

  template_dir = 'impsolv'
  def __init__(self, d):
    self.param_dict = d
    d['template_dir'] = 'impsolv'
    self.stages = self.param_dict.get('stages')
    self.inp_fn_list=[ 'toppar.str'
                     , 'step1_pdbreader.inp'
                     , 'step2_implicit.inp'
                     , '1d5r.pdb'
                     , 'pten_charged_membrane_start.inp'
                     , 'pten_charged_membrane_1.inp'
                     , 'pten_non_phos_gbsw_start.inp'
                     , 'pten_non_phos_gbsw_1.inp'
                     , 'pten_phos_start.inp'
                     , 'pten_phos_1.inp'
                     ]
    self.default_out_ext = '.pdb'
    s0in = d.get('charmm_s0in') or 'step1_pdbreader.inp'
    s1in = d.get('charmm_s1in') or 'step2_implicit.inp'
 

    self.stage_list = [[{'cmd': 'charmm-medium', 'si': s0in, 'a': 'toppar.str'}],
                       [{'cmd': 'charmm-medium', 'si': s1in, 'a': 'toppar.str'}],
                      ]                     

  @staticmethod
  def preprocess(d, local_paths):
    
    d['dest_dir']        = os.path.join(local_paths['resource_prefix'],
                                        d['dbname'])


    session_dir = d['session_dir']    
    #pdb_path = os.path.join('/damsl/projects/molecules/data/PTEN/data/josh/output', seq_id, seq_id+'.B99990001.pdb')
   
    output_prefix = os.path.join(local_paths['resource_prefix'], local_paths['io_dir'], session_dir, d['run_dir'])
    
    conn = mddb_utils.get_dbconn(d['dbname'],
                                 d['dbuser'],
                                 d['dbhost'],
                                 d['dbpass'])

    if not os.path.exists(output_prefix):
      os.makedirs(output_prefix)
    
    print 'output_prefix: ', output_prefix
    
    if ('sim_input_id' in d):
      cur = conn.cursor()
      cur.execute('select sim_input_fn from sim_inputs where sim_input_id = {0}').format(d['sim_input_id'])

      sim_input_fn = cur.fetchone()['sim_input_fn']
      d['charmm_s0in'] = sim_input_fn + '_start,in'
      d['charmm_s1in'] = sim_input_fn + '_1,in'


    if ('pdb_id' not in d):   
      # need to do structure prediction. First insert sequence into db
      insert_query = "SELECT sequence_insert(%s)"
      aaseq = d['aaseq']
      fasta = ""
      
      cur1 = conn.cursor()
      cur1.execute(insert_query,(fasta,)) 
      rows = cur1.fetchall()
      
       
       
      seqid = str(rows[0][0])
      d['seqid'] = seqid 
      cur1.close()
      conn.commit()
      
      fasta = ">" + seqid + "\n" + aaseq
      update_query = "UPDATE sequences_josh SET fasta=%s where seqid=%s"
      cur2 = conn.cursor()
      cur2.execute(update_query,(fasta, int(seqid)))
      cur2.close()
      conn.commit() 
            

      #  create modeller pdb
      # create file names for later use 
      fastaFileName = output_prefix + seqid + ".fasta"
      f = open(fastaFileName,'w+')
      f.write(fasta)
      f.close()     
      pirFileName = output_prefix + seqid +'.ali'
      # TODO check this
      pdbFileName = os.path.join(local_paths['template_prefix'].format("/home/webapp",d['user']),ImpSolvPipeline.template_dir,'1d5r.pdb')
      alignFileName = output_prefix + seqid +'-p60484'+ '.ali'
      # convert fasta file to pir file
      MT.fastaToPir(fastaFileName,pirFileName)
      # align pir to pdb file
      MT.alignPirToPdb(pirFileName,pdbFileName,seqid,alignFileName)
      # insert alignment to DB
      MT.dbInsertAlignment(conn, seqid, alignFileName)
      # build homology models
      MT.doAutoModel(alignFileName,pdbFileName,seqid)
      # move output to proper directory
      pdbfile = seqid + ".B99990001.pdb"
      f = open(pdbfile,'r')
      pdb = f.read()
      f.close()
               
      insert_query = "SELECT pdb_insert(%s,%s,%s)"
      
      cur1 = conn.cursor()
      cur1.execute(insert_query,(seqid,pdb,"modeller"))
      rows = cur1.fetchall()
      id = rows[0][0]
      cur1.close()
      conn.commit()
      d['pdb_id'] = id
       

    pdb_id = d['pdb_id']
    cur = conn.cursor()
    cur.execute('select pdb from seq_pdb where pdb_id = {0}'.format(pdb_id));
    pdb_content = cur.fetchone()[0]

 
    with open(os.path.join(output_prefix, 'protein0.pdb'), 'w') as ofp:
      ofp.write(pdb_content)

    cmd0 = 'python /damsl/projects/molecules/software/tools/pdb2pqr/pdb2pqr-1.8/pdb2pqr.py ' +\
          '--ffout=charmm --ff=charmm {0} {1}'.format('protein0.pdb', 'protein1.pdb')  
    # preprocess pdb file
    cmd1 = 'python /damsl/projects/molecules/data/PTEN/Full/Charmm/ACE/charmm-gui-model1/add_segid_PTEN.py ' +\
           '-f protein1.pdb -o protein2.pdb -p pdb -s ' +\
           '/damsl/projects/molecules/data/PTEN/Full/Charmm/ACE/charmm-gui-model1/resids.dat' 
    
    cmd2 = "sed 's/ ATOM/ATOM/' protein2.pdb | " +\
           "sed 's/CD1 ILE/CD  ILE/' | " +\
           "sed 's/1CBDISU/ CB CYS/' | " +\
           "sed 's/1SGDISU/ SG CYS/' > protein.pdb "   
 
    print cmd0
    subprocess.call(cmd0, shell=True, cwd = output_prefix)
    print cmd1
    subprocess.call(cmd1, shell=True, cwd = output_prefix)
    print cmd2
    subprocess.call(cmd2, shell=True, cwd = output_prefix)    

    #os.remove(os.path.join(output_prefix, 'protein0.pdb'))
    #os.remove(os.path.join(output_prefix, 'protein1.pdb'))
    #os.remove(os.path.join(output_prefix, 'protein2.pdb'))

  def run_substage(self, output_prefix, substage):

    in_fn = substage.get('i') or substage.get('si')

    out_fns = substage.get('o')

    cmd   = self.param_dict[substage.get('cmd')]

    log_fn = os.path.splitext(in_fn)[0] + '.log'
    subdir = substage.get('d') or ""
    workdir = os.path.join(output_prefix,subdir)

    print '==============================================================='
    print 'in_fn:', in_fn
    print 'output_prefix:', output_prefix
    print 'subdir:', subdir
    print 'workdir:', workdir
    print 'log_fn:', log_fn


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


    time.sleep(0.1)
    subprocess.call("tail -n 60 {0}".format(log_fn), shell=True, cwd = output_prefix,
                      stderr=subprocess.STDOUT)

    print out_fns
    print map(lambda out_fn: os.path.isfile(os.path.join(workdir, out_fn)), out_fns)

    #return all(map(lambda out_fn: os.path.isfile(os.path.join(workdir, out_fn)), out_fns))    
    return True

def test():
  d = {'resource_prefix': '/damsl/projects',
       'io_dir':          'resource_io',
       'seq_id':          '3', 
       'session_dir':     'dummy_session2',
       'run_dir':         'run_1',
       'template_prefix': '/home/chao/june27/damsl.mddb/templates/',
       'template_dir':    'impsolv',  
       'charmm':          '/damsl/projects/molecules/software/MD/Charmm/c36b1/exec/gnu/charmm',
      }

  local_paths = resources.LocalResource.get_paths()

  ImpSolvPipeline.preprocess(d, local_paths)

  pipeline = ImpSolvPipeline(d)
  pipeline.run()

def test2():
  d = {'resource_prefix': '/damsl/projects',
       'io_dir':          'resource_io',
       'session_dir':     'dummy_session3',
       'run_dir':         'run_1',
       'template_prefix': '/home/jbw/mddbnew2/damsl.mddb/templates/',
       'template_dir':    'impsolv',
       'dbname':          'mddb_bdslss3',
       'dbhost':          'mddb',
       'dbuser':          'jbw',
       'dbpass':          '/home/jbw/pass.out',
       'aaseq':            'MTAIIKEIVSRNKRRYQEDGFDLDLTYIYPNIIAMGFPAERLEGVYRNNIDDVVRFLDSKHKNHYKIYNLCAERHYDTAKFNCRVAQYPFEDHNPPQLELIKPFCEDLDQWLSEDDNHVAAIHCKAGKGRTGVMICAYLLHRGKFLKAQEALDFYGEVRTRDKKGVTIPSQRRYVYYYSYLLKNHLDYRPVALLFHKMMFETIPMFSGGTCNPQFVVCQLKVKIYSSNSGPTRREDKFMYFEFPQPLPVCGDIKVEFFHKQNKMLKKDKMFRFWVNTFFIPGPEETSEKVENGSLCDQEIDSICSIERADNDKEYLVLTLTKNDLDKANKDKANRYFSPNFKVKLYFTKTVEEPSNPEASSSTSVTPDVSDNEPDHYRYSDTTDSDPENEPFDEDQHTQITKV'
      }

  try:
    import resources
    local_paths = resources.LocalResource.get_paths()
  
    ImpSolvPipeline.preprocess(d, local_paths)
  
    pipeline = ImpSolvPipeline(d)
    pipeline.run()
  except Exception as e:
     print e

if __name__ == '__main__':
  test2()

