import functools, itertools, os, os.path, shutil, string, sys, re
from optparse import OptionParser
from itertools import izip
import cStringIO

def compose_csv(job_id, results):
  output = cStringIO.StringIO()
  map(lambda x: output.write("{0},{1}\n".format(job_id, ','.join(map(str,x)))), results)
  print output.getvalue()
  output.seek(0)
  return output


def parse_fred(fred_out_fn):
  with open(fred_out_fn, "r") as ifs:
    ifs.readline()
    return map(lambda l: l.split()[1:], ifs.readlines())
  return None

def parse_fred2(fred_out_fn):
  results = []
  with open(fred_out_fn, "r") as ifs:
    for line in ifs:
      if ':' in line:
        continue
      results.append(line.split()[1:])
  return results


def parse_dock6(dock6_out_fn):
  results = []
  with open(dock6_out_fn, 'r') as ifs:
   gs = None
   gvdw = None
   ges = None
   inen = None
   for line in ifs:
      if 'Grid Score:' in line:
        gs = float(line.split(':')[1])
      if 'Grid_vdw:' in line:
        gvdw = float(line.split(':')[1])
      if 'Grid_es:' in line:
        ges = float(line.split(':')[1])
      if 'Int_energy:' in line:
        inen = float(line.split(':')[1])
        if gs != None and gvdw != None and ges != None and inen != None:
          results.append([gs,gvdw,ges,inen])
          gs = None
          gvdw = None
          ges = None
          inen = None
  return results

def parse_dock6_amber(amber_out_fn):
  results = []
  with open(amber_out_fn, 'r') as ifs:
    bond = None  
    angle = None   
    dihedral = None   
    enb14 = None   
    eel14 = None  
    enb = None  
    eel = None  
    egb = None  
    econs = None   
    esurf = None    
    total = None
    amber_score = None
    complex = None
    receptor = None
    ligand = None
    for line in ifs:
      if  'bond:' in line:
        bond = float(line.split(':')[1])
      if  'angle:' in line:
        angle = float(line.split(':')[1])
      if  'dihedral:' in line:
        dihedral = float(line.split(':')[1])
      if  'enb14:' in line:
        enb14 = float(line.split(':')[1])
      if  'eel14:' in line:
        eel14 = float(line.split(':')[1])
      if  'enb:' in line:
        enb = float(line.split(':')[1])
      if  'eel:' in line:
        eel = float(line.split(':')[1])
      if  'egb:' in line:
        egb = float(line.split(':')[1])
      if  'econs:' in line:
        econs = float(line.split(':')[1])
      if  'esurf:' in line:
        esurf = float(line.split(':')[1])
      if  'Total:' in line:
        total = float(line.split(':')[1])
      if 'Amber Score:' in line:
        amber_score = float(line.split(':')[1])
      if 'complex:' in line:
        complex = float(line.split(':')[1])
      if 'receptor:' in line:
        receptor = float(line.split(':')[1])
      if 'ligand:' in line:
        ligand = float(line.split(':')[1])
        results.append([bond,angle,dihedral,enb14,eel14,enb,eel,
                        egb,econs,esurf,total,amber_score,complex,
                        receptor,ligand])
  return results

