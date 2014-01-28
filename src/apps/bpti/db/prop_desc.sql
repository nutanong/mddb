drop table if exists HydrogenBonds;

create table HydrogenBonds as
select distinct on (t,res_id2) * from (
  select distinct on (t,res_id1) * from (
    select A1.trj_id, A1.t,
           P1.residue_id as res_id1, 
           P2.residue_id as res_id2, 
           vec_length(A1.x-A2.x, A1.y-A2.y, A1.z-A2.z) as dist
    from trj_17705_ap A1, trj_17705_ap A2, ProteinAtoms P1, ProteinAtoms P2
    where P1.protein_id = 10 and P2.protein_id = 10 and
          A1.atom_id = P1.atom_id and A2.atom_id = P2.atom_id and
          P1.atom_name = 'HN' and P2.atom_name = 'O' and
          not P1.residue_id between P2.residue_id-1 and P2.residue_id+1 and
          A1.trj_id = A2.trj_id and A1.t = A2.t
  ) S1
  order by t,res_id1,dist
) S2
order by t,res_id2,dist;

drop table if exists SaltBridges;

create table SaltBridges as
select distinct on (t,res_id2) * from (
  select distinct on (t,res_id1) * from (
    select A1.trj_id, A1.t,
           P1.residue_id as res_id1, P1.atom_name as atom_name1,
           P2.residue_id as res_id2, P2.atom_name as atom_name2,
           vec_length(A1.x-A2.x, A1.y-A2.y, A1.z-A2.z) as dist
    from trj_17705_ap A1, ProteinAtoms P1, trj_17705_ap A2, ProteinAtoms P2
    where P1.protein_id = 10 and P2.protein_id = 10 and
          A1.atom_id = P1.atom_id and A2.atom_id = P2.atom_id and
          not P1.residue_id between P2.residue_id-1 and P2.residue_id+1 and
          ((P1.residue_name = 'LYS' and P1.atom_name = 'NZ') or 
           (P1.residue_name = 'ARG' and (P1.atom_name = 'NH1' or P1.atom_name = 'NH2'))) and
          (P2.residue_name = 'ASP' or P2.residue_name = 'GLU') and
          (P2.atom_name = 'OE1' or P2.atom_name = 'OE2') and
          A1.trj_id = A2.trj_id and A1.t = A2.t 
  ) S1
  order by t,res_id1,dist
) S2
order by t,res_id2,dist;
   

