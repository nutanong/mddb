import functools, itertools, os, os.path, shutil, string, sys, re
from optparse import OptionParser
from itertools import izip
import MDAnalysis, numpy.linalg
import traceback
import cStringIO
import math
import time

num_timesteps = pow(10,5)

# Enable debug output.
debug = False

# Default simulation package is CHARMM ('charmm').
# Alternative options are Amber ('amber') and Gromacs ('gromacs')
sim = "charmm"

# The different sections in the CHARMM PSF and PARAM files.
sections_charmm = ["atom", "bond", "angle", "di", "imp"]

# The different sections in the Amber PSF and PARAM files.
sections_amber_psf = ["atom",            "atom_name",            "atom_type",
                      "bond",            "bonds_inc_H",          "bonds_without_H",
                      "angle",           "angle_inc_H",          "angle_without_H",
                      "di",              "di_inc_H",             "di_without_H",
                      "imp",             "unique",               "uniqueimp",
                      "residue_label",   "residue_pointer",      "residue_temp",
                      "atom_type_index", "nonbonded_parm_index", "number_excluded_atoms"]

sections_amber_param = ["charge", "atom",
                        "bond", "psf_bond",
                        "angle", "psf_angle",
                        "di", "psf_dihedral",
                        "imp", "improper", "imph",
                        "bond_force_constant", "bond_equil_value",
                        "angle_force_constant", "angle_equil_value",
                        "dihedral_force_constant", "dihedral_periodicity", "dihedral_phase",
                        "lennard_jones_acoef", "lennard_jones_bcoef", "nb","phi","theta","nbond","psf_nb"]

# The different sections in the gromacs top and itp files.
sections_gromacs = ["atom", "bond", "angle", "di", "imp", "nb_partial", "pair", "atomtypes", "pairtypes"]


############################
#
# Constants

charmm_param_name_map = {'atom'  : ['atom'],
                         'bond'  : ['bond'],
                         'angle' : ['angle',     'theta'],
                         'di'    : ['dihedral',  'phi'],
                         'imp'   : ['improper',  'imph'],
                         'nb'    : ['nonbonded', 'nbond'] }

amber_file_name_map = {
    "atom"                    : "atom",
    "atom_name"               : "atom_name",
    "atom_type"               : "amber_atom_type",
    "atom_type_index"         : "atom_type_index",
    "nonbonded_parm_index"    : "nonbonded_parm_index",
    "charge"                  : "charge",
    "residue_label"           : "residue_label",
    "residue_pointer"         : "residue_pointer",
    "residue_temp"            : "residue_temp",
    "bonds_inc_H"             : "bonds_inc_hydrogen",
    "bonds_without_H"         : "bonds_without_hydrogen",
    "bond_force_constant"     : "bond_force_constant",
    "bond_equil_value"        : "bond_equil_value",
    "angle_inc_H"             : "angles_inc_hydrogen",
    "angle_without_H"         : "angles_without_hydrogen",
    "angle_force_constant"    : "angle_force_constant",
    "angle_equil_value"       : "angle_equil_value",
    "di_inc_H"                : "dihedrals_inc_hydrogen",
    "di_without_H"            : "dihedrals_without_hydrogen",
    "dihedral_force_constant" : "dihedral_force_constant",
    "dihedral_periodicity"    : "dihedral_periodicity",
    "dihedral_phase"          : "dihedral_phase",
    "lennard_jones_acoef"     : "lennard_jones_acoef",
    "lennard_jones_bcoef"     : "lennard_jones_bcoef",
    "bond"                    : "bond",
    "angle"                   : "angle",
    "di"                      : "dihedral",
    "nb"                      : "nonbonded",
    "number_excluded_atoms"   : "number_excluded_atoms",
    "improper"                : "improper",
    "imp"                     : "imp",
    "unique"                  : "unique",
    "uniqueimp"               : "uniqueimp",
    "psf_nb"                  : "psf_nb"
   }

gromacs_section_names = {
  "atom"        : "atom",
  "bond"        : "bond",
  "angle"       : "angle",
  "dihedral"    : "di",
  "improper"    : "imp", 
  "nonbond"     : "nb_partial",
  "pair"        : "pair",
  "atomtypes"   : "atomtypes",
  "pairtypes"   : "pairtypes",
  "atom_params" : "atom_params",
  "system"      : "system"
}

gromacs_improper_default_multiplicity = "0"

#############################
#
# Multi-file utilities

def zip_files(print_index, in_files, out_file, hdr):
  in_fs = map(lambda f: open(f, "r"), in_files)
  out_f = open(out_file, "w")
  add_output(out_f, hdr)
  for i, in_lines in enumerate(itertools.izip(*in_fs)):
    if i != 0:
      out = (str(i) if print_index else "")+(",".join([l.strip() for l in in_lines]))
      add_output(out_f, out+"\n")

def cat_files(in_files, out_file, hdr):
  out_f = open(out_file, "w")
  add_output(out_f, hdr)
  for f in in_files:
    x = open (f, "r")
    for index, line1 in enumerate(x):
       if index != 0:
         line = line1.strip()+"\n"
         add_output(out_f,line)


#####################
##                 ##
## DCD parsing     ############################################################
##                 ##
#####################

def write_trajectory_range(universe, start, stop, step, output_fn, trj_id):
  output_line = []
  count       = 0
  output_st   = cStringIO.StringIO()
  output_file = open(output_fn, "w") if output_fn <> None else None

  atomlist = list(universe.atoms)

  for ts in universe.trajectory[start:stop:step]:
    coords = universe.atoms.coordinates()
    timestep = ts.frame
    for atom in atomlist:
      atom_id = atom.number + 1
      x = coords[atom_id-1,0]
      y = coords[atom_id-1,1]
      z = coords[atom_id-1,2]
      if debug:
        print "Atom: {0} (t:{4}, x:{1}, y:{2}, z:{3})".format(atom_id, x, y, z, timestep)
      
      output_fields = map(lambda x: str(x), [trj_id,timestep, atom_id, x, y, z])
      if output_file == None:
        output_st.write("\t".join(output_fields)+"\n")

      else:
        output_line.append(",".join(output_fields)+"\n")
        count = count + 1
        if count%1000 == 0:
          output_file.writelines(output_line)
          output_file.flush()
          output_line=[]
  
  if output_file == None:
    output_st.seek(0)
  else:
    output_file.writelines(output_line)
    output_file.flush()
  
  return output_st if output_file == None else None;



def write_trajectory(universe, output_fn, trj_id):
  numAtoms = universe.atoms.numberOfAtoms() 

  print "Number of atoms in the PSF file is %d" % numAtoms
  print "Number of frames %d" % len(universe.trajectory)

  return write_trajectory_range(universe, 0, -1, 1, output_fn, trj_id)

# This method is called to extract the trajectory information from the DCD file.
def parse_trajectory_to_file(trj_id, top_fn, trj_fn, output_fn):
  universe = MDAnalysis.Universe(top_fn, trj_fn) 
  write_trajectory(universe, output_fn, trj_id)
 

def parse_trajectory_range_to_buffer(trj_id, top_fn, trj_fn, start, stop, step):
  print "Loading trajectory"
  universe = None
  try:
    universe = MDAnalysis.Universe(str(top_fn), str(trj_fn))
  except:
    traceback.print_stack()
    print sys.exc_info()[0]
    print top_fn, trj_fn

  return write_trajectory_range(universe, start, stop, step, None, trj_id) if universe <> None else None



def parse_trajectory_to_buffer(trj_id, top_fn, trj_fn):
  print "Loading trajectory"
  universe = None
  try:
    universe = MDAnalysis.Universe(str(top_fn), str(trj_fn)) 
  except:
    traceback.print_stack()
    print sys.exc_info()[0]
    print top_fn, trj_fn

  return write_trajectory(universe, None, trj_id) if universe <> None else None

def get_trjectory_length(top_fn, trj_fn):
  try:
    universe = MDAnalysis.Universe(str(top_fn), str(trj_fn))
  except:
    traceback.print_stack()
    print sys.exc_info()
    print top_fn, trj_fn

  return len(universe.trajectory)



def parse_trajectory_to_buffer_partial(trj_id, num_blocks, block_num, top_fn, trj_fn):
  print "Loading trajectory"
  universe = None
  start, stop, step = 0, -1, 1

  try:
    universe = MDAnalysis.Universe(top_fn, trj_fn)
    num_frames = len(universe.trajectory)
    frames_per_block = num_frames/num_blocks
    start = block_num * frames_per_block    
    stop = start + frames_per_block
    step = 10
    print "start, stop: ", start, stop
  except:
    print sys.exc_info()[0]
    print top_fn, trj_fn

  return write_trajectory_range(universe, start, stop, step, None, trj_id) if universe <> None else None


 
###########################
#
# Generic parsing helpers

output_lines = []

def is_partial_line(line):
  return all(map(lambda x: not(x in line), ["!", "*"])) \
            and line.strip().endswith("-")

def has_header(match_start, values_to_match, line):
  return any(map(
    lambda h: line.startswith(h) if match_start else (h in line),
    values_to_match))
  
def header(match_start, values_to_match, line):
  return filter(
    lambda h: line.startswith(h) if match_start else (h in line),
    values_to_match)

def buffer_output(line):
  global output_lines
  output_lines.append(line)

# Used to write output lines to the specified output line. output_lines is
# a global var which is reset after writing the lines to the output file.
def add_output(output_file, line):
  if line == '\n' or len(line) == 0:
    return

  global output_lines
  output_lines.append(line)
  output_file.writelines(output_lines)
  output_lines = []
  
# Used to reset output_lines global var.
def reset_output():
  global output_lines
  output_lines = []


###################
# Field handlers
# All field handlers should return a boolean indicating whether to search
# again for the section header.

# This method builds tuples of the given size from sequential fields
def tuple_of_seq_fields(id, tuple_size, fields, output_file):
  splits = []
  for i in range(tuple_size):
    splits.append(fields[i::tuple_size])
  tuples = zip(*splits)
  for t in tuples:
    if debug:
      print "{0}: {1}".format(id, ",".join(t))
    add_output(output_file, ",".join(t)+"\n")
  return True

# This method writes the fields obtained from line parser to respective output file.
def project_record_fields(id, projections, fields, output_file):
  out_fields = [fields[i] for i in projections]
  if debug:
    print "{0}: {1}".format(id, ",".join(out_fields))
  add_output(output_file, ",".join(out_fields)+"\n")
  return True

#added an extra field handler method for gromacs parsing
check_tuple = []

def project_record_fields_with_check(id, projections1, projections2, fields1, fields2,
                                     fields1_as_ids, id_type, output_file):
  global check_tuple
  out_fields = [id_type[fields1[i]] if fields1_as_ids else fields1[i] for i in projections1] \
                + [fields2[i] for i in projections2]

  if (id=='bond' or id=='angle'):
    curr_tuple     = map(lambda x:out_fields[x],range(0,len(projections1)))
    curr_tuple_rev = map(lambda x:out_fields[x],range(len(projections1)-1,-1,-1))
    if ( curr_tuple not in check_tuple and curr_tuple_rev not in check_tuple):
      check_tuple = check_tuple + [curr_tuple]
      if debug:
        print "{0}: {1}".format(id, ",".join(out_fields))
      add_output(output_file, ",".join(out_fields))

  elif (id=='dihedral' or id=='improper'):
    curr_tuple       = map(lambda x:out_fields[x],range(0,len(projections1)))
    curr_tuple_rev   = map(lambda x:out_fields[x],range(len(projections1)-1,-1,-1))
    curr_tuple_exch1 = [curr_tuple[0],curr_tuple[2],curr_tuple[1],curr_tuple[3]]
    curr_tuple_exch2 = [curr_tuple[3],curr_tuple[1],curr_tuple[2],curr_tuple[0]]
    if ( curr_tuple not in check_tuple and curr_tuple_rev not in check_tuple and
         curr_tuple_exch1 not in check_tuple and curr_tuple_exch2 not in check_tuple):
      check_tuple = check_tuple + [curr_tuple]
      if debug:
          print "{0}: {1}".format(id, ",".join(out_fields))
      add_output(output_file, ",".join(out_fields))

  return True


###############
# Line parsers.
#
# These accept an accumulated section string, a single line string, and
# a file handle on which to output data.
#
# All line parsers should return a boolean indicating whether to search
# again for the section header.

def direct_line_parser(fields_handler, section, line, output_file):
  return fields_handler(line, output_file), section
        
# This method parses a section like atom, bond, angle, dihedral etc.
# line_parser is a handle to functools.partial(parse_line, line_parser).
# hdr_parser is a handle to functools.partial(parse_line, hdr_parser) in which
# hdr_parser is the section header string like NATOM/NBOND/Bond/Atom etc.
def parse_section(input, output, output_hdr, match_start,
                  possible_vals, desired_val, hdr_parser, line_parser):
  partial_line = ""
  section = ""
  plist = []
  with open(input) as in_file:
    with open(output, "a") as output_file:
      next_line = in_file.next
      buffer_output(output_hdr+"\n") 
      some_found = found = False
      for line in in_file:
        if is_partial_line(line):
          partial_line += line
        else:
          line = partial_line + line
          partial_line = ""

          is_hdr = has_header(match_start, possible_vals, line)
          hdr_match = desired_val in header(match_start, possible_vals, line)
          
          if found and is_hdr and not(hdr_match):
            found = False

          elif not found and is_hdr and hdr_match:
            some_found = found = True
            try:
              if debug:
                print "Parsing header line {0}".format(line.strip())
              hdr_parser(line, output_file)
            except Exception:
              print "Failure parsing header, aborting section {0}".format(desired_val)
              break

          elif found:
            try:
              found, section = line_parser(section, line, output_file)
            except Exception:
              print "Failure parsing section {0} line, skipping: {1}".format(desired_val, line.strip())
              found = False

          elif debug:
            print "Skipping: \"{0}\"".format(line.strip())
      
      if not some_found:
        print "Warning:Failed to find {0} section in {1}".format(desired_val, input)
      elif debug:
        print "Found {0} section in {1}".format(desired_val, input)
      
      reset_output()

def section_to_csv(input, output, output_hdr, match_start, headers, id, line_parser):
  possible_hdrs, _ = zip(*headers.values())
  section_hdr      = headers[id][0]
  hdr_parser       = headers[id][1]
  parse_section(input, output, output_hdr, match_start,
                possible_hdrs, section_hdr, hdr_parser, line_parser)



########################
##                    ##
## CHARMM             #########################################################
##                    ##
########################

########################
#
# CHARMM line and field parsers

def charmm_section_line_parser(fields_handler, section, line, output_file):
  # sline contains the number of atoms/bonds/angles etc. mentioned in
  # the section header.
  sline = line.strip().partition("!")[0] 

  # Ignore empty and comment lines.
  if sline != "" and all(map(lambda x: not(sline.startswith(x)), ["!", "*",";"])):
      return fields_handler(sline.split(), output_file), section
  else:
      return True, section

def charmm_field_parser(fields_handler, id):
  return functools.partial(charmm_section_line_parser,
            functools.partial(fields_handler,id), "")  

##########################
#
# CHARMM section parsers

def charmm_section_to_csv(input, output, output_hdr, match_start,
                          headers, id, fields_handler):
  section_to_csv(input, output, output_hdr, match_start, headers, id,
    functools.partial(charmm_section_line_parser, fields_handler))

def charmm_tuple_section_to_csv(match_start, headers, id, tuple_size,
                                output_hdr, input, output):
  charmm_section_to_csv(input, output, output_hdr, match_start, headers, id,
    functools.partial(tuple_of_seq_fields, id, tuple_size))

def charmm_record_section_to_csv(match_start, headers, id, projections,
                                 output_hdr, input, output):
  charmm_section_to_csv(input, output, output_hdr, match_start, headers, id,
    functools.partial(project_record_fields, id, projections))


####################################
#
# PSF section parsing for CHARMM

def psf_header(id, fields, output_file): 
  if debug:
    print "CHARMM PSF #{0}: {1}".format(id, fields[0])

psf_hdr_parser = lambda x: charmm_field_parser(psf_header, x)

psf_headers = dict(
  # Creating key value pairs viz. atom = "!NATOM"
  map(lambda x: [x,("!N"+x.upper(), psf_hdr_parser(x))],
    ["atom", "bond", "theta", "phi", "imphi", "don", "acc", "nb"]))

# Section parsers
charmm_psf_parsers = {
 "record" : [ ("atom", [0,5,4,2,3,1], 
               "atom_id,atom_ty,atom_nm,res_id,res_nm,seg_nm") ],

 "tuple" : [ ("bond", 2, "atom_id1,atom_id2"),
             ("theta", 3, "atom_id1,atom_id2,atom_id3"),
             ("phi", 4, "atom_id1,atom_id2,atom_id3,atom_id4"),
             ("imphi", 4, "atom_id1,atom_id2,atom_id3,atom_id4") ]
}

def parse_charmm_psf(psf, outputs):
  for sec_type, sections in charmm_psf_parsers.iteritems():
    if sec_type == "record":
      for parse_args in sections:
        output = outputs[parse_args[0]] 
        args = parse_args + (psf, output)
        charmm_record_section_to_csv(False, psf_headers, *args)
    
    elif sec_type == "tuple":
      for parse_args in sections:
        output = outputs[parse_args[0]] 
        args = parse_args + (psf, output)
        charmm_tuple_section_to_csv(False, psf_headers, *args)
    
    else:
      print "Invalid parser type: {0}".format(sec_type)


#########################################
#
# Paramfile section parsing for CHARMM

def prm_header(id, fields, output_file):
  if debug:
    print "PRM {0}: {1}".format(id, fields[0])

prm_hdr_parser = lambda x: charmm_field_parser(prm_header, x)

prm_headers = dict(
  map(lambda x: [x, (x.upper(), prm_hdr_parser(x))],
      ["bond", "angle", "theta", "dihedral", "phi",
       "improper", "imph", "nonbonded", "nbond"]))

ang_args = (range(5),     "atom_ty1,atom_ty2,atom_ty3,Ktheta,theta0")
di_args  = (range(7),     "atom_ty1,atom_ty2,atom_ty3,atom_ty4,Kchi,n,delta")
imp_args = (range(7),     "atom_ty1,atom_ty2,atom_ty3,atom_ty4,Kchi,n,delta")
nb_args  = ([0, 2, 3],    "atom_ty,rmin/2,epsilon")

charmm_prm_parsers = {
  "record" : [("bond",         range(4), "atom_ty1,atom_ty2,Kb,b0"),
              ("angle",)     + ang_args,
              ("theta",)     + ang_args,
              ("dihedral",)  + di_args,
              ("phi",)       + di_args,
              ("improper",)  + imp_args,
              ("imph",)      + imp_args,
              ("nonbonded",) + nb_args,
              ("nbond",)     + nb_args],
  "tuple"  : []
}

def parse_charmm_prm(prm, outputs):
  for sec_type, sections in charmm_prm_parsers.iteritems():
    if sec_type == "record":
      for parse_args in sections:
        output = outputs[parse_args[0]] 
        args = parse_args + (prm, output)
        charmm_record_section_to_csv(True, prm_headers, *args)

    elif sec_type == "tuple":
      for parse_args in sections:
        output = outputs[parse_args[0]] 
        args = parse_args + (prm, output)
        charmm_tuple_seq_section_to_csv(True, prm_headers, *args)

    else:
      print "Invalid parser type: {0}".format(sec_type)


#####################################
#
# Topology file parsing for CHARMM

rtf_residue_id = ""
rtf_residue_atom_group = []
rtf_atoms = {}

def rtf_header(fields, output_file):
  global rtf_residue_id
  rtf_residue_id = fields[1]

rtf_headers = dict(map(lambda x: [x, (x.upper(), rtf_header)], ["resi"]))

def rtf_next_atom_group(output_file):
  global rtf_residue_atom_group
  if rtf_residue_atom_group != []:
    data = "\n".join(map(lambda x: ",".join([rtf_residue_id]+x), rtf_residue_atom_group))
    add_output(output_file,data+"\n")
  rtf_residue_atom_group = []

def rtf_residue_line(fields, output_file):
  global rtf_residue_atom_group
  # Handle residue atoms statefully, by tracking atoms per residue group
  # across multiple lines.
  # Other elements of the residue can be handled similarly here as needed.
  reset = False
  if "ATOM" in fields[0]:
    rtf_residue_atom_group.append([fields[2], fields[1], fields[3]])
  else:
    # Reset the search for residues as soon as we're finished with atom groups.
    if not("GROU" in fields[0]):
      reset = True
    rtf_next_atom_group(output_file)
  return not(reset)
  
def rtf_residues_to_csv(rtf, output):
  rtf_resfld_parser   = functools.partial(rtf_residue_line)
  rtf_hdr_matches, _  = zip(*rtf_headers.values())
  rtf_sec_hdr         = rtf_headers["resi"][0]
  parse_section(rtf, output, "res_id,atom_ty,atom_nm,charge", True,
                rtf_hdr_matches, rtf_sec_hdr,
                functools.partial(charmm_section_line_parser, rtf_header, ""),
                functools.partial(charmm_section_line_parser, rtf_resfld_parser))


def rtf_atom_line(fields, output_file):
  global rtf_atoms
  reset = False
  if "ATOM" in fields[0]:
    nm,ty,ch = fields[1:4]
    if (rtf_residue_id,ty,nm) not in rtf_atoms:
      rtf_atoms[rtf_residue_id,ty,nm] = ch
  else:
    # Reset the search for residues as soon as we're finished with atom groups.
    if not("GROU" in fields[0]):
      reset = True
  return not(reset)

def rtf_atoms_to_csv(rtf, output):
  global rtf_atoms
  rtf_atoms          = {}
  rtf_atfld_parser   = functools.partial(rtf_atom_line)
  rtf_hdr_matches, _ = zip(*rtf_headers.values())
  rtf_sec_hdr        = rtf_headers["resi"][0]
  output_hdr         = "residue_nm,atom_ty,atom_nm,charge"
  
  parse_section(rtf, output, output_hdr, True,
                rtf_hdr_matches, rtf_sec_hdr,
                functools.partial(charmm_section_line_parser, rtf_header, ""),
                functools.partial(charmm_section_line_parser, rtf_atfld_parser))
  
  if rtf_atoms != {}:
    with open(output, "a") as output_file:
      buffer_output(output_hdr+"\n")
      data = "\n".join(map(lambda ((rn,ty,nm),ch): ",".join([rn,ty,nm,ch]), rtf_atoms.items()))
      add_output(output_file, data+"\n")
  rtf_atoms = {}



########################
##                    ##
## AMBER              #########################################################
##                    ##
########################

######################
#
# Amber line parsers

# Amber parsing globals

section_counter = 0
section_counts  = {}

num_atoms            = ""
num_atom_types       = ""
num_bonds_inc_H      = ""
num_bonds_without_H  = ""
num_angles_inc_H     = ""
num_angles_without_H = ""
atom_type_list       = []
nbparm_index_list    = []
atom_type_index_list = []

impoutput     = None
uniqueoutput  = None
uniqimpoutput = None


def amber_output_dihedral_tuple(tuple, output_file):
  abs_int_list = lambda l,mask: \
    ','.join([str(abs(int(x))) if y else x for x,y in zip(l,mask)])+"\n"

  mask4        = [True, True, True, True]

  tuple_file   = None
  tuple_str    = abs_int_list(tuple[0:4], mask4)
  elem3, elem4 = int (tuple[2]), int (tuple[3])
  
  if elem4<0:
    tuple_file = open(impoutput, "a")
  elif elem3<0:
    tuple_file = open(uniqueoutput, "a")
  else: 
    tuple_file = output_file

  if tuple_file <> None:
    add_output(tuple_file, tuple_str)

  if elem4<0 and elem3<0:
    amber_uniqueimp_ids_file = open(uniqimpoutput,"a")
    add_output(amber_uniqueimp_ids_file, tuple_str)

def amber_tuple_of_seq_fields(tuple_id, tuple_size, fields, output_file):

  abs_int_list = lambda l,mask: \
    ','.join([str((abs(int(x))/3)+1) if y else x for x,y in zip(l,mask)])+"\n"

  mask2  = [True, True]
  mask3  = [True, True, True]
  mask21 = [True, True, False]
  mask31 = [True, True, True, False]

  if tuple_id == "Bond_force_constant"     or tuple_id == "Bond_equil_value" or \
     tuple_id == "Angle_force_constant"    or tuple_id == "Angle_equil_value" or \
     tuple_id == "Dihedral_force_constant" or \
     tuple_id == "Dihedral_periodicity"    or tuple_id == "Dihedral_phase" or \
     tuple_id == "Lennard_jones_acoef"     or tuple_id == "Lennard_jones_bcoef" or \
     tuple_id == "Residue_label" or tuple_id == "Residue_pointer":
    fields = fields.split()
    for t in fields:
      add_output(output_file, t+"\n")

  elif tuple_id == "Atom_name":
    fields = fields.split()
    for t in fields:
      atom4 = len(t)/4
      if len(t) > 4:
        for i in range(1,atom4+1):
          tuple = t[4*(i-1):(4*i)]
          add_output(output_file, tuple+"\n")
        add_output(output_file, t[4*atom4:]+"\n")
      else:
        add_output(output_file, t+"\n")

  elif tuple_id == "Atom_type_index":
    fields = fields.split()
    for i in fields:
      atom_type_index_list.append(i)
      add_output(output_file, i+"\n")

  elif tuple_id == "Amber_atom_type":
    fields = fields.split()
    for i in fields:
      atom_type_list.append(i)
      add_output(output_file, i+"\n") 
        
  elif tuple_id == "Charge":
    fields = fields.split()
    for t in fields:
      chargematch = re.match("^([-]*[0-9]+[.][0-9]+)E([+-][0-9]+).*", t)
      if chargematch:
        mantissa = chargematch.group(1)
        exponent = chargematch.group(2)
        t = (float(mantissa))*(pow(10, int(exponent)))
        t = t/18.2223
        add_output(output_file, str(t)+"\n")

  elif tuple_id == "Nonbonded_parm_index":
    fields = fields.split() 
    for i in fields:
      nbparm_index_list.append(i)
      add_output(output_file, i+"\n")

  elif tuple_id == "Bonds_inc_hydrogen" or tuple_id == "Bonds_without_hydrogen":
    fields =  fields.split()
    count = len(fields)
    for i in range(0,count,3):
      add_output(output_file, abs_int_list(fields[i:i+2], mask2))

  elif tuple_id == "Angles_inc_hydrogen" or tuple_id == "Angles_without_hydrogen":
    fields =  fields.split()
    count = len(fields)
    for i in range(0,count,4):
      add_output(output_file, abs_int_list(fields[i:i+3], mask3))

  elif tuple_id == "Dihedrals_inc_hydrogen" or tuple_id == "Dihedrals_without_hydrogen":
    fields = fields.split()
    tuple1 = fields[:5]
    amber_output_dihedral_tuple(tuple1, output_file)

    if len(fields) >= 10:
      tuple2 = fields[5:]
      amber_output_dihedral_tuple(tuple2, output_file)

  return True



###############################################
# Helper functions for Amber param file parsing

# This method finds out the relevant section by matching every line with required number of hyphens and other criteria.
def has_amber_section(line, section_name, section_id):
  global section_counter
  tuples = line.split('-')

  tsplit = lambda x: sum([tuples[i].split() for i in range(0,x)], [])
  valid  = lambda l,s,pos: \
    len(line) > l and all(map(lambda x: (line[x[0]] == '-') == x[1], pos)) and \
    section_name.lower() == s and section_counter == section_id

  # Matching bond section.  
  if valid(5, "bond", [(2,True), (5,False)]):
    line = tsplit(2)[:4]     
  # Matching angles/theta section.  
  elif valid(8, "angle", [(2,True), (5,True), (8,False)]): 
    line = tsplit(3)[:5] 
  # Matching dihedrals/phi section. Need changes to correctly parse dihedrals in the end of the file. #AskYanif/Tom  
  elif valid(10, "di", [(2,True), (5,True), (8,True)]):
    line = tsplit(4)[:7] 
  # Matching impropers section.
  elif valid(10, "imp", [(2,True), (5,True), (8,True)]):
    line = tsplit(4)[:7] 
  # Matching nonbonded section.  
  elif re.match("\s{2}[A-Z].[ ]{10}[0-9][.][0-9]+.*",line) and section_name.lower() == "nb":
    nonbondedtuples = line.split()
    line = nonbondedtuples[:3]
  else:
    line = ""
  return line


##########################
#
# Amber section parsers
def amber_append_section(desired_val, section, line):
  new_section = section
  if "FORMAT" in line:
    return new_section

  if  desired_val == "%FLAG BONDS_INC_HYDROGEN"      or desired_val == "%FLAG BONDS_WITHOUT_HYDROGEN" or \
      desired_val == "%FLAG ANGLES_WITHOUT_HYDROGEN" or desired_val == "%FLAG ANGLES_INC_HYDROGEN"    or \
      desired_val == "%FLAG ATOM_TYPE_INDEX"         or desired_val == "%FLAG NONBONDED_PARM_INDEX" :
    new_section += line
  return new_section

def amber_line_parser(desired_val, fields_handler, section, line, output_file):
  new_section = amber_append_section(desired_val, section, line)
  # Ignore empty and comment lines.
  if not(line.startswith("%FORMAT")):
      return fields_handler(line, output_file), new_section
  else:
      return True, new_section

def amber_whole_section_parser(desired_val, count, section_parser, section, line, output_file):
  new_section = amber_append_section(desired_val, section, line)
  found = True
  if len(new_section.split()) == count:
    found = section_parser(new_section, output_file)
    new_section = ""
  return found, new_section

# This method parses each section in the Amber Param file by reading each line and
# then passing the record to record line parser to write the sections to files.
def amber_parse_section(match_start, id, count, input, output, output_hdr, line_parser):
  global section_counter
  partial_line = ""
  with open(input) as in_file:
    with open(output, "a") as output_file:
      next_line = in_file.next
      buffer_output(output_hdr+"\n") # Writes output header to the file.
      for line in in_file:
        if line.strip() == "": # If a blank file is found it indicates a section will start.
          line = next_line() # Skip the blank line.
          section_counter += 1

        # Checking the new line if it is a section. After the method call the
        # tuples corresponding to some section will be returned in line.
        line = has_amber_section(line, id, count) 
        if line != '':
          if debug:
            print "Parsing section {0} line {1}".format(section_counter, line)
          line_parser("", line, output_file) # Write the tuples to the file.
      reset_output()


def amber_section_to_csv(match_start, id, count, input, output, output_hdr, fields_handler):
  amber_parse_section(match_start, id, count, input, output, output_hdr,
                      functools.partial(direct_line_parser, fields_handler))

def amber_record_section_to_csv(match_start, id, count, projections, output_hdr, input, output):
  amber_section_to_csv(match_start, id, count, input, output, output_hdr, \
                       functools.partial(project_record_fields, id.capitalize(), projections))
                        # The last argument is a partial function in which the first 2 arguments of
                        # project_record_fields are being set to id.capitalize() and projections.
                        # Note that the partial function is just being defined by setting the arguments
                        # and not being called.

def amber_tuple_section_to_csv(match_start, headers, id, tuple_size, output_hdr, count, input, output):
  field_parser = functools.partial(amber_tuple_of_seq_fields, id.capitalize(), tuple_size)
  if count > 0:
    line_parser = functools.partial(amber_whole_section_parser, headers[id][0], count, field_parser)
  else:
    line_parser = functools.partial(amber_line_parser, headers[id][0], field_parser)
  section_to_csv(input, output, output_hdr, match_start, headers, id, line_parser)


#########################
#
# Amber param parsing

amber_param_parsers = {
  "record"  : [ ("bond",  1, range(4),         "atom_ty1,atom_ty2,Kb,b0"),
                ("angle", 2, range(5),         "atom_ty1,atom_ty2,atom_ty3,Ktheta,theta0"),
                ("di",    3, range(4)+[6,4,5], "atom_ty1,atom_ty2,atom_ty3,atom_ty4,Kchi,n,delta"),
                ("imp",   4, range(4)+[6,4,5], "atom_ty1,atom_ty2,atom_ty3,atom_ty4,Kchi,n,delta"),
                ("nb",    7, range(3),         "atom_ty,rmin/2,epsilon")]
}

def parse_amber_param(prm, outputs):
  global section_counter
  for sec_type, sections in amber_param_parsers.iteritems():
    if sec_type == "record":
      for parse_args in sections:
        output = outputs[parse_args[0]]
        args = parse_args + (prm, output)
        amber_record_section_to_csv(False, *args)
        section_counter = 0


#################################
#
# PSF section parsing for Amber

def psf_amber_header(id, fields, output_file):
  if debug:
    print "AMBER PSF #{0}: {1}".format(id, fields[0])

psf_amber_headers = dict(map( \
  lambda x: [x,("%FLAG "+x.upper(),functools.partial(psf_amber_header,x))],
  ["charge", "mass",
   "tree_chain_classification",
   "atom_name", "amber_atom_type", "atom_type_index",
   "nonbonded_parm_index",
   "bonds_inc_hydrogen",      "bonds_without_hydrogen",
   "angles_inc_hydrogen",     "angles_without_hydrogen",
   "dihedrals_inc_hydrogen",  "dihedrals_without_hydrogen",
   "excluded_atoms_list",
   "bond_force_constant",     "bond_equil_value",
   "angle_force_constant",    "angle_equil_value",
   "dihedral_force_constant", "dihedral_periodicity", "dihedral_phase",
   "scee_scale_factor",
   "lennard_jones_acoef",     "lennard_jones_bcoef",
   "residue_label",           "residue_pointer", "number_excluded_atoms"]))

# Section parsers
amber_psf_parsers = {
  "tuple"  : [ ("atom_name",                  1, "atom_ty"),
               ("amber_atom_type",            1, "atom_nm"),
               ("atom_type_index",            1, "atom_type_index"),
               ("nonbonded_parm_index",       1, "nonbonded_parm_index"),
               ("charge",                     1, "charge"),
               ("residue_label",              1, "residue_label"),
               ("residue_pointer",            1, "residue_pointer"),
               ("bonds_inc_hydrogen",         3, "atom_id1,atom_id2,index"),
               ("bonds_without_hydrogen",     3, "atom_id1,atom_id2,index"),
               ("bond_force_constant",        1, "bond_force_constant"),
               ("bond_equil_value",           1, "bond_equil_value"),
               ("angles_inc_hydrogen",        4, "atom_id1,atom_id2,atom_id3,index"),
               ("angles_without_hydrogen",    4, "atom_id1,atom_id2,atom_id3,index"),
               ("angle_force_constant",       1, "angle_force_constant"),
               ("angle_equil_value",          1, "angle_equil_value"),
               ("dihedrals_inc_hydrogen",     5, "atom_id1,atom_id2,atom_id3,atom_id4,index"),
               ("dihedrals_without_hydrogen", 5, "atom_id1,atom_id2,atom_id3,atom_id4,index"),
               ("dihedral_force_constant",    1, "dihedral_force_constant"),
               ("dihedral_periodicity",       1, "dihedral_periodicity"),
               ("dihedral_phase",             1, "dihedral_phase"),
               ("lennard_jones_acoef",        1, "lennard_jones_acoef"),
               ("lennard_jones_bcoef",        1, "lennard_jones_bcoef") ] 
}

def psf_amber_pointers(psf):
  global section_counts
  global atom_type_index_list
  global atom_type_list
  global nbparm_index_list
  global num_atom_types
  global num_atoms
  global num_bonds_inc_H
  global num_bonds_without_H
  global num_angles_without_H
  global num_angles_inc_H
  with open(psf) as in_file:
    for line in in_file:
      if re.match("^%FLAG POINTERS.*", line):
        line = in_file.next()
        line = in_file.next()
        line = line.split()
        atom_type_index_list = []
        atom_type_list       = []
        nbparm_index_list    = []
        num_atoms            = line[0]
        num_atom_types       = line[1]
        num_bonds_inc_H      = line[2]
        num_bonds_without_H  = line[3]
        num_angles_inc_H     = line[4]
        num_angles_without_H = line[5]
        section_counts       = {
          "bonds_inc_hydrogen"      : int(num_bonds_inc_H) * 3,
          "bonds_without_hydrogen"  : int(num_bonds_without_H) * 3,
          "angles_inc_hydrogen"     : int(num_angles_inc_H) * 4,
          "angles_without_hydrogen" : int(num_angles_without_H) * 4,
          "atom_type_index"         : int(num_atoms),
          "nonbonded_parm_index"    : int(num_atom_types) * int(num_atom_types)
        }

        print "Number of atoms in PSF file: %s"                           % line[0]
        print "Number of distinct atom types in PSF file: %s"             % line[1]
        print "Number of distinct bond ids with H in PSF file: %s"        % line[2]
        print "Number of distinct bond ids without H in PSF file: %s"     % line[3]
        print "Number of distinct angles ids with H in PSF file: %s"      % line[4]
        print "Number of distinct angles ids without H in PSF file: %s"   % line[5]
        print "Number of distinct dihedral ids with H in PSF file: %s"    % line[6]
        print "Number of distinct dihedral ids without H in PSF file: %s" % line[7]

def parse_amber_psf(psf, outputs):
  for sec_type, sections in amber_psf_parsers.iteritems():
    if sec_type == "tuple":
      for parse_args in sections:
        if debug:
          print "Parsing Amber PSF {0}".format(parse_args[0])
        count = 0
        if parse_args[0] in section_counts:
          count = section_counts[parse_args[0]]
        output = outputs[parse_args[0]]
        args = parse_args + (count, psf, output)
        amber_tuple_section_to_csv(True, psf_amber_headers, *args)

      if debug:
        print "Finished parsing Amber PSF tuple sections"

    else:
      print "Invalid parser type: {0}".format(sec_type)

amber_psf_id_csv_deps = {
  "bond"  : ("cat", False,
             "atom_id1,atom_id2\n",
             ["bonds_inc_hydrogen", "bonds_without_hydrogen"]),
  
  "angle" : ("cat", False,
             "atom_id1,atom_id2,atom_id3\n",
             ["angles_inc_hydrogen", "angles_without_hydrogen"]),

  "di"    : ("cat", False, 
             "atom_id1,atom_id2,atom_id3,atom_id4\n",
             ["dihedrals_inc_hydrogen", "dihedrals_without_hydrogen"])
}

amber_psf_params_csv_deps = {
  "psf_bond"     : ("zip", False,
                    "bond_force_constant,bond_equil_value\n",
                    ["bond_force_constant", "bond_equil_value"]),
  
  "psf_angle"    : ("zip", False,
                    "angle_force_constant,angle_equil_value\n",
                    ["angle_force_constant", "angle_equil_value"]),
  
  "psf_dihedral" : ("zip", False,
                    "dihedral_force_constant,dihedral_periodicity,dihedral_phase \n",
                    ["dihedral_force_constant", "dihedral_periodicity", "dihedral_phase"])
}


def psf_amber_write_csvs(amber_files, mapped_files, csv_dict):
  for k, (out_type, print_index, hdr, sec_keys) in csv_dict.iteritems():
    in_files = map(lambda x: mapped_files[x], sec_keys)
    if out_type == "zip":
      zip_files(print_index, in_files, amber_files[k], hdr)
    elif out_type == "cat":
      cat_files(in_files, amber_files[k], hdr)
    else:
      print "Invalid output type: {0}".format(out_type)

def psf_amber_id_param_csv(psf_files, param_files, mapped_files):
  psf_amber_write_csvs(psf_files, mapped_files, amber_psf_id_csv_deps)
  psf_amber_write_csvs(param_files, mapped_files, amber_psf_params_csv_deps)


def psf_amber_atom_fields_to_csv(add_counter, header, in_fn, out_fn):
  files = [open(fn,"r") for fn in in_fn]
  output_f = open(out_fn, "w")
  add_output(output_f, header)
  for i, lines in enumerate(izip(*files)):
    if i != 0:
      fields = ([str(i)] if add_counter else [])+[l.strip() for l in lines] + ["blah"]
      line = ",".join(fields)+"\n"
      add_output(output_f, line)

def psf_amber_atom_params_to_csv(psf_file, param_file, mapped_files):
  fn     = [mapped_files[k] for k in \
            ["atom_name", "amber_atom_type", "charge", "residue_temp"]]
  header = "atom_ty,atom_nm,charge,res_id,res_nm,seg_nm\n"
  psf_amber_atom_fields_to_csv(False, header, fn, param_file["atom"])

def psf_amber_atom_ids_to_csv(psf_file, param_file, mapped_files):
  fn     = [mapped_files[k] for k in \
            ["atom_name", "amber_atom_type", "residue_temp"]]
  header = "atom_id,atom_ty,atom_nm,res_id,res_nm\n"
  psf_amber_atom_fields_to_csv(True, header, fn, psf_file["atom"])


def psf_amber_atom_residue_ids_to_csv(amber_files):
  in_fs   = map(lambda x: open(amber_files[x], "r"), ["residue_label", "residue_pointer"])
  out_f   = open(amber_files["residue_temp"], "w")
  header  = "res_id,res_nm\n"
  add_output(out_f, header)  
  in_lines   = map(lambda f: [l.strip() for l in f.readlines()], in_fs)
  min_length = min(len(in_lines[0]), len(in_lines[1]))
  for i in range(1, min_length):
    l = int (in_lines[1][i])
    u = int (in_lines[1][i+1]) if i != (min_length-1) else int (num_atoms) + 1
    for j in range(l,u):
      line = (str(i+1) + "," + str(in_lines[0][i]) + "\n")
      add_output(out_f,line)

def psf_amber_nonbonded_params_to_csv(amber_files):
  in_fs  = map(lambda x : open(amber_files[x], "r"), ["lennard_jones_acoef", "lennard_jones_bcoef"])
  out_f  = open(amber_files["psf_nb"], "w")
  header = "atom_id1,atom_id2,atom_ty1,atom_ty2,acoeff,bcoeff\n"
  add_output(out_f, header)  
  num        = len(atom_type_index_list)
  num_types  = int (num_atom_types) 
  in_lines   = map(lambda f: [l.strip() for l in f.readlines()][1:], in_fs)

  print "Amber PSF NBP length {0}".format(str(num))

  #for i in range(0, num_types):
  #  for j in range(0, num_types):
  #    index = nbparm_index_list[num_types * i + j]
  #    a = in_lines[0][int(index) - 1]
  #    b = in_lines[1][int(index) - 1]
  #    line = str(i+1)+","+str(j+1)+","+str(atom_type_list[i])+","+str(atom_type_list[j])+","+str(a) +","+str(b) +"\n"
  #    add_output(out_f,line)

  for i in range(0, num):
    for j in range(0,num):
      value = ((num_types*(int(atom_type_index_list[i]) -1)) + int (atom_type_index_list[j]))
      if value <= (num_types * num_types):
        index = nbparm_index_list[value - 1]   
        a = in_lines[0][int(index) - 1]
        b = in_lines[1][int(index) - 1]
        line = str(i+1)+","+str(j+1)+","+str(atom_type_list[i])+","+str(atom_type_list[j])+","+str(a) +","+str(b) +"\n"
        add_output(out_f,line)


######################################################
# 
# PDB file parsing for Amber to get residue information

def rtf_amber_header(fields, output_file):
  global rtf_residue_id
  rtf_residue_id = fields[0]

rtf_amber_headers = dict(map(lambda x: [x, (x.upper(), rtf_header)], ["int"]))

def rtf_amber_atom_line(fields, output_file):
  global rtf_atoms
  reset = False
  if len(fields) > 1:
    if not("CORR" in fields[0]) and len(fields) > 10:
      (ty,nm),ch = fields[1:3], fields[10]
      if (rtf_residue_id,ty,nm) not in rtf_atoms:
        rtf_atoms[rtf_residue_id,ty,nm] = ch
  
  elif "DONE" in fields[0]:
      reset = True
  
  return not(reset)

def rtf_amber_atoms_to_csv(rtf, output):
  global rtf_atoms
  rtf_atoms          = {}
  rtf_atfld_parser   = functools.partial(rtf_amber_atom_line)
  rtf_hdr_matches, _ = zip(*rtf_amber_headers.values())
  rtf_sec_hdr        = rtf_amber_headers["int"][0]
  output_hdr         = "residue_nm,atom_ty,atom_nm,charge"
    
  parse_section(rtf, output, output_hdr, False,
                rtf_hdr_matches, rtf_sec_hdr,
                functools.partial(charmm_section_line_parser, rtf_amber_header, ""),
                functools.partial(charmm_section_line_parser, rtf_atfld_parser))
  
  if rtf_atoms != {}:
    with open(output, "a") as output_file:
      buffer_output(output_hdr+"\n")
      data = "\n".join(map(lambda ((rn,ty,nm),ch): ",".join([rn,ty,nm,ch]), rtf_atoms.items()))
      add_output(output_file, data+"\n")
  rtf_atoms = {}


########################
##                    ##
## GROMACS            #########################################################
##                    ##
########################

def strip_gromacs_index(args):
  new_args = []
  for id, proj, hdr in args:
    hdr_fields = hdr.split(",")
    new_args += \
      [(id, proj[:-1], ",".join(hdr_fields[:-1]))] \
      if "index_variable" in hdr_fields or "funct" in hdr_fields else \
      [(id, proj, hdr)]

  return new_args

#######################################
#
# Topology file parsing for Gromacs

dihedral_sec_counter = 0

def gromacs_top_header(id, fields, output_file):
  if debug:
    print "TOP #{0}: {1}".format(id, fields[0])

def gromacs_improper_top_header(id, fields, output_file):
  global dihedral_sec_counter
  dihedral_sec_counter += 1
  if debug:
    print "TOP #{0}: {1}".format(id, fields[0])

def gromacs_improper_top_record(id, projections, fields, output_file):
  if dihedral_sec_counter <= 1:
    return False
  return project_record_fields(id, projections, fields, output_file)

def gromacs_impropers_top_section(match_start, headers, id, projections,
                                  output_hdr, input, output):
  section_to_csv(input, output, output_hdr, match_start, headers, id,
                 functools.partial(charmm_section_line_parser, 
                   functools.partial(gromacs_improper_top_record, id, projections)))

hdr_parser = lambda x: \
  (x, charmm_field_parser(gromacs_top_header, x)) \
  if x <> "improper" else \
  ("dihedral", charmm_field_parser(gromacs_improper_top_header,x))

gromacs_top_headers = dict(
    map(lambda x: [x,hdr_parser(x)],
    ["atom", "bond", "pair", "angle", "dihedral","improper","system"]))

# Section parsers
gromacs_top_parsers = {
 "record" : [ ("atom",        [0,1,4,2,3], "atom_id,atom_ty,atom_nm,res_id,res_nm"),
              ("atom_params", [1,4,3,6],   "atom_ty,atom_nm,res_nm,charge"),
              ("bond",        [0,1,3],     "atom_id1,atom_id2,index_variable"),
              ("angle",       [0,1,2,4],   "atom_id1,atom_id2,atom_id3,index_variable"),
              ("dihedral",    [0,1,2,3,5], "atom_id1,atom_id2,atom_id3,atom_id4,index_variable"),
              ("improper",    [0,1,2,3,5], "atom_id1,atom_id2,atom_id3,atom_id4,index_variable"),
              ("pair",        [0,1,2],     "atom_id1,atom_id2,funct")]
}

def gromacs_parse_top(input, outputs, strip_index):
  parsers = {'record': strip_gromacs_index(gromacs_top_parsers["record"])} \
            if strip_index else gromacs_top_parsers

  for sec_type, sections in parsers.iteritems():
    if sec_type == "record":
      for parse_args in sections:
        output = outputs[gromacs_section_names[parse_args[0]]]
        args   = parse_args + (input, output) \
                 if parse_args[0] <> "atom_params" else \
                 ("atom", parse_args[1], parse_args[2]) + (input, output)
        if parse_args[0] == "improper":
          gromacs_impropers_top_section(False, gromacs_top_headers, *args)
        else:
          charmm_record_section_to_csv(False, gromacs_top_headers, *args)

    else:
      print "Invalid parser type: {0}".format(sec_type)

#for splitting the pair_ids file into nonbond_params_ids and pairtypes_ids file
def split_pair_ids(inputs,outputs):
  di_pair=[]
  with open(inputs[1]) as in_file:
    for line in in_file:
      if line.find("atom_id")==-1:
        tokens=line.split(',')
        di_pair=di_pair+[[tokens[0],tokens[3]]]


  with open(inputs[0]) as in_file:
    with open(outputs[0], "a") as output_file_nb:
      with open(outputs[1],"a") as output_file_pt:
        for line in in_file:
          if line.find("atom_id")==-1:
            tokens=line.split(',')
            if([tokens[0],tokens[1]] in di_pair or [tokens[1],tokens[0]] in di_pair ):
              output_file_pt.writelines(line)
            else:
              output_file_nb.writelines(line)
          else:
            output_file_nb.writelines(line)
            output_file_pt.writelines(line)


#############################################
#
# Gromacs RTP (residue topology) parser

def gromacs_rtp_header(id, fields, output_file):
  if debug:
    print "RTP #{0}: {1}".format(id, fields[0])

header_func = lambda x: charmm_field_parser(gromacs_rtp_header,x)
  
gromacs_rtp_sections        = ["atom", "bond", "angle", "dihedral","improper","exclusion"]
gromacs_rtp_section_matches = gromacs_rtp_sections + [ "bondedtypes" ]
gromacs_rtp_headers         = dict(map(lambda x: [x,("[ "+x+"s ]", header_func(x))], gromacs_rtp_sections))

# Specialized methods for atom parsing
def gromacs_atom_rtp_header(id, fields, output_file):
  global rtf_residue_id
  rtf_residue_id = fields[1]
  if debug:
    print "RTP #{0}: {1}".format(id, " ".join(fields))

def gromacs_atom_rtp_hdr_parser(fields_handler, section, line, output_file):
  global rtf_residue_id
  sline = line.strip().partition("!")[0] 
  if sline != "" and all(map(lambda x: not(sline.startswith(x)), ["!", "*",";"])):
    if sline.startswith("[") \
        and all(map(lambda x: not(x in sline), gromacs_rtp_section_matches)):
      fields_handler(sline.split(), output_file)
    else:
      rtf_residue_id = ""

def gromacs_atom_rtp_record(id, projections, fields, output_file):
  global rtf_residue_id
  
  is_section = lambda f: \
    all(map(lambda x: not(x in f), gromacs_rtp_section_matches))

  if rtf_residue_id <> "":
    # Handle sequential headers
    if "[" in fields:
      r = (" ".join(fields) == "[ atoms ]")
      return r
    else:
      # Print if we have a valid residue, and this is not a header.
      out_fields = [fields[i] for i in projections]
      out_fields = [rtf_residue_id] + out_fields
      if debug:
        print "{0}: {1}".format(id, ",".join(out_fields))
      add_output(output_file, ",".join(out_fields)+"\n")
      return True

  elif "[" in fields and is_section(fields[1]):
    rtf_residue_id = fields[1]
    return True
  
  # All other cases are invalid atom sections
  return False

def gromacs_atom_rtp_section(match_start, headers, id, projections,
                              output_hdr, input, output):
  section_to_csv(input, output, output_hdr, match_start, headers, id,
                 functools.partial(charmm_section_line_parser, 
                   functools.partial(gromacs_atom_rtp_record, id, projections)))

atom_header_func = lambda x: \
  functools.partial(gromacs_atom_rtp_hdr_parser, functools.partial(gromacs_atom_rtp_header,x), "")

#gromacs_atom_rtp_matches = [("residue", "[")] + map(lambda x: (x," [ "+x+"s ]"), gromacs_rtp_sections)
gromacs_atom_rtp_matches = [("atom", "[")]
gromacs_atom_rtp_headers = dict(map(lambda (x,y): [x,(y, atom_header_func(x))], gromacs_atom_rtp_matches))

gromacs_rtp_parsers = {
 "record" : [ ("atom", [1,0,2],"residue_nm,atom_ty,atom_nm,charge"),
              ("bond",[0,1,2],"atom_id1,atom_id2,index_variable"),
              ("angle",[0,1,2,3],"atom_id1,atom_id2,atom_id3,index_variable"),
              ("dihedral",[0,1,2,3,4],"atom_id1,atom_id2,atom_id3,atom_id4,index_variable"),
              ("improper",[0,1,2,3,4],"atom_id1,atom_id2,atom_id3,atom_id4,index_variable")]
}

def gromacs_parse_rtp(input, outputs, strip_index):
  parsers = {'record': strip_gromacs_index(gromacs_rtp_parsers["record"])} \
            if strip_index else gromacs_rtp_parsers

  for sec_type, sections in parsers.iteritems():
    if sec_type == "record":
      for parse_args in sections:
        output = outputs[gromacs_section_names[parse_args[0]]]
        args   = parse_args + (input, output)
        if parse_args[0] == "atom":
          gromacs_atom_rtp_section(False, gromacs_atom_rtp_headers, *args)
        else:
          charmm_record_section_to_csv(False, gromacs_rtp_headers, *args)

    else:
      print "Invalid parser type: {0}".format(sec_type)


#############################################
#
# ITP file(parameter) parsing for Gromacs

def gromacs_prm_header(id, fields, output_file):
  if debug:
    print "PRM {0}: {1}".format(id, fields[0])

prm_header = lambda x: charmm_field_parser(gromacs_prm_header,x)
gromacs_prm_headers = \
  {'bond'    :  ('Bond type code',                    prm_header("bond")),
   'angle'   :  ('Bond-angle type code',              prm_header("angle")),
   'dihedral':  ('Dihedral-angle type code',          prm_header("dihedral")),
   'improper':  ('Improper dihedral-angle type code', prm_header("improper")),
   'atomtypes': ('atomtypes',                         prm_header("atomtypes")),
   'nonbond':   ('nonbond',                           prm_header("nonbond")),
   'pairtypes': ('pairtypes',                         prm_header("pairtypes"))
   }

#Section parsers
gromacs_prm_parsers = {
  "record" : [("bond",     [1,3,2],   "index_func,Kb,b0"),
              ("angle",    [1,3,2],   "index_func,Ktheta,theta0"),
              ("dihedral", [1,2,3,4], "index_func,Kchi,n,delta"),
              ("improper", [1,2,3],   "index_func,Kchi,n,delta")]
}

gromacs_nbprm_parsers = {
  "record" : [("atomtypes",[0,5,6],   "atom_ty,C6,C12"),
              ("nonbond",  [0,1,3,4], "atom_ty,atom_ty2,C6,C12"),
              ("pairtypes",[0,1,3,4], "atom_ty,atom_ty2,C6,C12")]
}


def gromacs_improper_prm_record(id, projections, fields, output_file):
  out_fields = [fields[i] for i in projections]
  # Inject value for multiplicity
  out_fields.insert(len(out_fields)-1, gromacs_improper_default_multiplicity) 
  if debug:
    print "{0}: {1}".format(id, ",".join(out_fields))
  add_output(output_file, ",".join(out_fields)+"\n")
  return True

def gromacs_impropers_prm_section(match_start, headers, id, projections,
                              output_hdr, input, output):
  section_to_csv(input, output, output_hdr, match_start, headers, id,
                 functools.partial(charmm_section_line_parser, 
                   functools.partial(gromacs_improper_prm_record, id, projections)))


def gromacs_parse_prm(gromacs_parsers, input, outputs):
  for sec_type, sections in gromacs_parsers.iteritems():
    if sec_type == "record":
      for parse_args in sections:
        output = outputs[gromacs_section_names[parse_args[0]]]
        args   = parse_args + (input, output)
        if parse_args[0] == "improper":
          gromacs_impropers_prm_section(False, gromacs_prm_headers, *args)
        else:
          charmm_record_section_to_csv(False, gromacs_prm_headers, *args)

    else:
      print "Invalid parser type: {0}".format(sec_type)


gromacs_final_prm_parsers = {
  "record" : [("bond",     [0,1],[1,2],       "atom_ty1,atom_ty2,Kb,b0"),
              ("angle",    [0,1,2],[1,2],     "atom_ty1,atom_ty2,atom_ty3,Ktheta,theta0"),
              ("dihedral", [0,1,2,3],[1,2,3], "atom_ty1,atom_ty2,atom_ty3,atom_ty4,Kchi,n,delta"),
              ("improper", [0,1,2,3],[1,2,3], "atom_ty1,atom_ty2,atom_ty3,atom_ty4,Kchi,n,delta")]
}

#combine the id file and inter-param file to get the final param file
def combine_param_files(id,projection1,projection2,header_fields,input1_as_ids,id_type,input1,input2,output):
  global check_tuple
  check_tuple = []
  
  mk_args = lambda line1,line2: \
    (id, projection1, projection2, line1.split(','), line2.split(','), \
     input1_as_ids, id_type, output_file)
  
  with open(output, "a") as output_file:
    output_file.writelines(header_fields + "\n")
    with open(input1) as in_file1:
      for line1 in in_file1:
        index_func1 = line1.split(',')[len(projection1)].strip()
        if (index_func1.find("index") == -1):  #skip the header line
          with open(input2) as in_file2:
            for line2 in in_file2:
              index_func2 = line2.split(',')[0].strip()
              if (index_func1 == index_func2):
                project_record_fields_with_check(*(mk_args(line1,line2)))
                break


def combine_nb_and_pt_param_files(inputs,output):
  with open(output,"a") as out_file:
    with open(inputs[0]) as in_file1:
      for line in in_file1:
        out_file.writelines(line)
    with open(inputs[1]) as in_file2:
      for line in in_file2:
        if (line.find("atom")==-1):
          out_file.writelines(line)


def make_param_files(input1_as_ids, inputs1, inputs2, outputs):
  id_type = {}
  if input1_as_ids:
    with open(inputs1["atom"]) as atom_id:
      for line in atom_id:
        if (line.split(',')[0].strip() != "atom_id"):
          id_type[line.split(',')[0].strip()]=line.split(',')[1].strip()

  shutil.copy(inputs1["atom"], outputs[gromacs_section_names["atom"]])

  for sec_type, sections in gromacs_final_prm_parsers.iteritems():
    if sec_type == "record":
      for parse_args in sections:
        input1 = inputs1[gromacs_section_names[parse_args[0]]]
        input2 = inputs2[gromacs_section_names[parse_args[0]]]
        output = outputs[gromacs_section_names[parse_args[0]]]
        args   = parse_args + (input1_as_ids, id_type, input1, input2, output)
        combine_param_files(*args)


def make_nbparam_files(outputs, pairwise):
  c6_idx  = 2 if pairwise else 1
  c12_idx = 3 if pairwise else 2
  header  = ",".join((["atom_ty1", "atom_ty2"] if pairwise else ["atom_ty"]) + ["sigma" , "epsilon"])+"\n"
  for src in outputs :
    with open(src) as in_file:
      with open(outputs[src], "a") as output_file:
        output_file.writelines(header)
        for line in in_file:
          if line.find("atom") == -1 and line != "\n":
            tokens=line.split(',')
            if float(tokens[c6_idx]) == 0 or float(tokens[c12_idx]) == 0:
              epsilon = sigma = 0
            else:
              sigma   = pow((float(tokens[c6_idx])/float(tokens[c12_idx])),(1.0/6))
              epsilon = float(tokens[c12_idx])/(4*pow(sigma,6))
            atom_ty = ([tokens[0], tokens[1]] if pairwise else [tokens[0]])+[str(sigma),str(epsilon)]
            output_file.writelines(",".join(atom_ty)+"\n")



#######################
##                   ##
## Toplevel          ##########################################################
##                   ##
#######################

# The program starts here. We add the options and the arguments parsing functionality here.
#

def validate_sim_mode(sim, mode):
  if sim.lower() <> "amber" and sim.lower() <> "charmm" and sim.lower() <> "gromacs":
    parser.error("Invalid simulation packkage: {0}".format(sim)) 

  if not (mode.lower() == "param" or mode.lower() == "protein" or mode.lower() == "trj"):
    parser.error("Invalid parse mode: {0}".format(mode)) 

def run_parser(sim, mode, 
               psf, dcd, param, nbparam, topology,
               in_dir, out_dir, out_prefix, trj_id):
  # Input files.
  # Assigning full path names for psf, dcd, prm and rtf variables.
  (psf, dcd, prm, rtf) = \
    map(lambda x: os.path.join(in_dir,x), [psf, dcd, param, topology])

  # Extract simulated system from the CHARMM trajectory.
  if sim == "charmm":
    
    # outputs list contains all the different output file names.
    outputs = \
      map(lambda x: os.path.join(out_dir, out_prefix+"_"+x+".csv"),
          (["charmm_"+"pos"] + 
           map(lambda x: "charmm_"+x+"_ids", sections_charmm) +  
           map(lambda x: "charmm_"+x+"_params", sections_charmm+["nb"]))) 

    charmm_positions = outputs[0]
    param_offset     = 1+len(sections_charmm)
    id_outputs       = dict(zip(sections_charmm, outputs[1:param_offset]))
    param_outputs    = dict(zip(sections_charmm+["nb"], outputs[param_offset:]))

    charmm_param_outputs = \
      dict(list(itertools.chain( \
        *[map(lambda n: (n,v), charmm_param_name_map[k]) \
          for k,v in param_outputs.iteritems()])))

    if mode.lower() == "trj":
      print "Calling CHARMM trajectory parsing function"
      parse_trajectory_to_file(trj_id, psf, dcd, charmm_positions)      

    elif mode.lower() == "protein":
      print "Calling CHARMM PSF parsing functions"
      id_name_map = {'atom' : 'atom', 'bond' : 'bond', 'angle' : 'theta',
                     'di' : 'phi', 'imp' : 'imphi' }

      charmm_id_outputs = dict([(id_name_map[k], v) for k,v in id_outputs.iteritems()])
      parse_charmm_psf(psf, charmm_id_outputs)

    else:
      # Extract constants from the CHARMM parameter file.
      print "Calling CHARMM Param parsing functions"  
      parse_charmm_prm(prm, charmm_param_outputs)

      # Extract charges from the CHARMM topology file.
      print "Calling CHARMM topology parsing functions"
      charmm_atom_params = outputs[6]
      rtf_atoms_to_csv(rtf, charmm_atom_params)

  # Extract simulated system from the Amber trajectory.
  elif sim == "amber":

    outputs_amber = \
      map(lambda x: os.path.join(out_dir, out_prefix+"_"+x+".csv"),
          (["amber_"+"pos"] + 
           map(lambda x: "amber_"+x+"_params", sections_amber_param+["nb"]) +  
           map(lambda x: "amber_"+x+"_ids",    sections_amber_psf))) 
    
    amber_positions = outputs_amber[0]
    psf_offset      = 2+len(sections_amber_param)
    param_outputs   = dict(zip(sections_amber_param, outputs_amber[1:psf_offset]))
    psf_outputs     = dict(zip(sections_amber_psf, outputs_amber[psf_offset:]))

    outputs1, outputs2 = (psf_outputs, param_outputs) if mode.lower() == "protein" \
                         else (param_outputs, psf_outputs)
    
    amber_outputs      = dict( [(amber_file_name_map[k], outputs1[k]) \
                                  if k in outputs1 else \
                                (amber_file_name_map[k], outputs2[k]) \
                                for k in amber_file_name_map.keys()])

    if mode.lower() == "trj":
      print "Calling Amber trajectory parsing functions"
      parse_trajectory_to_file(trj_id, psf, dcd, amber_positions)    

    elif mode.lower() == "protein": 
      global impoutput
      global uniqueoutput
      global uniqimpoutput

      impoutput     = amber_outputs["imp"]
      uniqueoutput  = amber_outputs["unique"]
      uniqimpoutput = amber_outputs["uniqueimp"]

      print "Calling Amber PSF/global file parsing functions"
      psf_amber_pointers(psf)
      parse_amber_psf(psf, amber_outputs)
      psf_amber_atom_residue_ids_to_csv(amber_outputs)
      psf_amber_atom_ids_to_csv(psf_outputs, param_outputs, amber_outputs)
      psf_amber_id_param_csv(psf_outputs, param_outputs, amber_outputs)
      #psf_amber_nonbonded_params_to_csv(amber_outputs)
      #psf_amber_atom_params_to_csv(psf_outputs, param_outputs, amber_outputs)

    else:
      # Extract constants from the Amber trajectory.
      print "Calling Amber param parsing functions"
      parse_amber_param(prm,param_outputs)

      print "Calling Amber topology parsing functions"
      amber_atom_params = amber_outputs["atom"]
      rtf_amber_atoms_to_csv(rtf, amber_atom_params)

    
  ###### Extract simulated system from the Gromacs trajectory.
  elif sim == "gromacs":

    # Assigning full path names of PSF, DCD, PARAM and TOPOLOGY files to gro,trr, prm, nbprm and top variables.
    (gro, trr, prm, nbprm, top) = \
      map(lambda x: os.path.join(in_dir,x), [psf, dcd, param, nbparam, topology])

    # outputs list contains all the different output file names.
    outputs = map(lambda x: os.path.join(out_dir, out_prefix+"_"+x+".csv"),
            (["gromacs_pos"] + 
             map(lambda x: "gromacs_"+x+"_ids",       sections_gromacs) +
             map(lambda x: "gromacs_"+x+"_param_tmp", sections_gromacs) +
             map(lambda x: "gromacs_"+x+"_params",    sections_gromacs) +
             map(lambda x: "gromacs_"+x+"_global",    sections_gromacs)))

    combined_nb_param_output  = os.path.join(out_dir,out_prefix+"_gromacs_nb_params.csv") 
   
    ns                        = len(sections_gromacs)
    off                       = zip(range(1, 3*ns+2, ns),range(ns+1, 4*ns+2, ns))
    param_offset              = 1+ns
    gromacs_positions         = outputs[0]
    
    dicts                     = [ dict(zip(sections_gromacs, outputs[off[i][0]:off[i][1]])) for i in range(0,4) ]
    id_outputs                = dicts[0]
    id_outputs["atom_params"] = outputs[19]
    inter_param_outputs       = dicts[1]
    param_outputs             = dicts[2]
    global_outputs            = dicts[3]
    
    if mode.lower() == "trj":
      print "Calling gromacs trajectory parsing function"
      parse_trajectory_to_file(trj_id, gro, trr, gromacs_positions)
    
    elif mode.lower() == "protein":
      print "Calling gromacs TOP parsing functions"
      gromacs_parse_top(top, id_outputs, True)
      nb_id_inputs  = [id_outputs[gromacs_section_names["pair"]],
                       id_outputs[gromacs_section_names["dihedral"]]]
      nb_id_outputs = [id_outputs[gromacs_section_names["nonbond"]],
                       id_outputs[gromacs_section_names["pairtypes"]]]
      split_pair_ids(nb_id_inputs,nb_id_outputs)

    else:
      print "Calling gromacs RTP parsing functions"
      gromacs_parse_rtp(top, global_outputs, False)

      # Extract constants from the gromacs parameter file.
      print "Calling gromacs Param parsing functions"
      gromacs_parse_prm(gromacs_prm_parsers, prm, inter_param_outputs)
      make_param_files(False, global_outputs, inter_param_outputs, param_outputs)

      print "Calling gromacs Nonbonded Param parsing functions"
      gromacs_parse_prm(gromacs_nbprm_parsers, nbprm, inter_param_outputs)

      #for conversion from C6 & C12 to epsilon & sigma
      nb_inkeys, nb_pairwise = ["atomtypes"], False # ["nonbond", "pairtypes"], True
      nb_keys                = map(lambda x: gromacs_section_names[x], nb_inkeys)
      nb_inter_param_inputs  = map(lambda x: inter_param_outputs[x], nb_keys)
      nb_param_outputs       = map(lambda x: param_outputs[x], nb_keys)

      make_nbparam_files(dict(zip(nb_inter_param_inputs, nb_param_outputs)), nb_pairwise)
      if nb_pairwise:
        combine_nb_and_pt_param_files(nb_param_outputs, combined_nb_param_output)
      else:
        shutil.copy(nb_param_outputs[0], combined_nb_param_output)
        
  else:
    print "Please enter either CHARMM or AMBER or GROMACS in the --sim or -s option"



if __name__ == '__main__':
  usage = "Usage: %prog [options] <output file prefix> <trj id>"
  parser = OptionParser(usage=usage)

  parser.add_option("-s", "--sim", type="string", dest="sim",
                    default="charmm",
                    help="select the simulator format from CHARMM/Amber/Gromacs", metavar="#SIM")
  
  parser.add_option("-m", "--mode", type="string", dest="mode",
                    default="trj",
                    help="set the parsing mode [param|protein|trj]", metavar="#MODE")

  parser.add_option("-i", "--inputdir", type="string", dest="inputdir",
                    default="alanine_dipeptide/",
                    help="set input directory", metavar="#INDIR")
    
  parser.add_option("-o", "--outputdir", type="string", dest="outputdir",
                    default="csv/",
                    help="set output directory", metavar="#OUTDIR")

  parser.add_option("-j", "--trjtop", type="string", dest="trjtop",
                    default="coor/alad_ace_dims_al_gen_250.psf",
                    help="set trajectory topology file", metavar="#TRJTOP")

  parser.add_option("-c", "--trj", type="string", dest="trj",
                    default="trj/tst_this_1.dcd",
                    help="set trajectory file", metavar="#TRJ")

  parser.add_option("-r", "--param", type="string", dest="param",
                    default="inp/param19-1.2.inp",
                    help="set param file", metavar="#PARAM")
  
  parser.add_option("-n", "--nbparam", type="string", dest="nbparam",   #added an extra option for gromacs nonbonded param file
                    default=" ",
                    help="set nbparam file", metavar="#NBPARAM")
  
  parser.add_option("-t", "--topology", type="string", dest="topology",
                    default="inp/toph19.rtf",
                    help="set topology file", metavar="#TOP")

  (options, args) = parser.parse_args()

  if len(options.sim) < 1:
    parser.error("No simulation package specified")
  
  if len(args) < 1:
    parser.error("no output file prefix specified")

  if len(args) < 2:
    parser.error("no trajectory id specified")

  sim     = options.sim.lower()
  mode    = options.mode

  validate_sim_mode(sim, mode)

  in_dir     = options.inputdir
  out_dir    = options.outputdir
  out_prefix = args[0]
  trj_id     = args[1]
  
  if not os.path.exists(in_dir):
    parser.error("invalid input directory "+in_dir)
    
  if not os.path.exists(out_dir):
    os.mkdir(out_dir)
  
  run_parser(sim, mode,
             options.trjtop, options.trj,
             options.param, options.nbparam, options.topology, 
             in_dir, out_dir, out_prefix, trj_id)
