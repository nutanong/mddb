
drop function if exists reinit_config_tables();
create or replace function reinit_config_tables()
returns bigint as $$
begin
  truncate table ConfigVDWNonBonded;
  truncate table ConfigDihedrals;
  truncate table ConfigAngles;
  truncate table ConfigBonds;
  truncate table ConfigAtoms;
  --truncate table AtomTypeIDs;
end
$$ language plpgsql;


drop function if exists reinit_protein_tables();
create or replace function reinit_protein_tables()
returns bigint as $$
begin
  truncate table ProteinElecNonBonded;
  truncate table ProteinFeatures;
  truncate table ProteinDihedrals;
  truncate table ProteinAngles;
  truncate table ProteinBonds;
  truncate table ProteinAtoms;
end
$$ language plpgsql;

drop function if exists insert_new_protein(text, text, text, text);
create or replace function insert_new_protein(
  new_protein_cname  text,
  new_protein_name   text, 
  new_protein_pdb_id text,
  new_protein_fn     text)
returns bigint as $$
declare
  new_protein_id bigint;
begin
  lock table Proteins in exclusive mode;
  insert into Proteins (protein_cname, protein_name, protein_pdb_id, protein_fn) 
    values (new_protein_cname, new_protein_name, new_protein_pdb_id, new_protein_fn);
  select max(protein_id) into new_protein_id from Proteins;

  return new_protein_id;
end
$$ language plpgsql;

create or replace function insert_new_config(
  sim_val       text,
  ff_val        text)
returns bigint as $$
declare
  new_config_id bigint;
begin
  lock table SimConfigs in exclusive mode;
  insert into SimConfigs (simulator, ff)
    values (sim_val, ff_val);
  select max(config_id) into new_config_id from SimConfigs;

  return new_config_id;
end
$$ language plpgsql;


----------------------------------------------------------------------------------------------
-- Config Loader
----------------------------------------------------------------------------------------------


drop function if exists load_config(bigint, text);
create or replace function load_config(
  c_id bigint,
  file_prefix text
)
returns void as $$
begin

  perform read_config_csv_files(file_prefix);

  perform load_config_atoms(c_id);
  perform load_config_bonds(c_id);
  perform load_config_angles(c_id);
  perform load_config_dihedrals(c_id);
  perform load_config_nonbondedpairs(c_id);

  --delete from VDW_NonBonded NB where NB.rmin = 0;
end
$$language plpgsql;

drop function if exists read_config_csv_files(text);
create or replace function read_config_csv_files(
  file_prefix text)
returns void as $$
begin

  drop table if exists AtomParams;

  create temp table AtomParams (
    residue_name  text,
    atom_type     text,
    atom_name     text,
    charge        float
  );

  execute 'copy AtomParams (residue_name, atom_type, atom_name, charge) ' ||
          ' from ' || quote_literal(file_prefix || '_atom_params.csv') ||
          ' with csv header';

  drop table if exists BondParams;
  create temp table BondParams (
    atom_ty1     text,
    atom_ty2     text,
    bond_const   float,
    bond_length  float
  );  

  execute 'copy BondParams ' ||
          ' from ' || quote_literal(file_prefix || '_bond_params.csv') ||
          ' with csv header';

  drop table if exists AngleParams;
  create temp table AngleParams (
    atom_ty1     text,
    atom_ty2     text,
    atom_ty3     text,
    angle_const  float,
    angle        float
  );

  execute 'copy AngleParams ' ||
          ' from ' || quote_literal(file_prefix || '_angle_params.csv') ||
          ' with csv header';

  drop table if exists DihedralParams;
  create temp table DihedralParams (
    atom_ty1     text,
    atom_ty2     text,
    atom_ty3     text,
    atom_ty4     text,
    force_const  float,
    n            float,
    delta        float 
  );

  execute 'copy DihedralParams ' ||
          ' from ' || quote_literal(file_prefix || '_di_params.csv') ||
          ' with csv header';


  drop table if exists ImproperParams;
  create temp table ImproperParams (
    atom_ty1     text,
    atom_ty2     text,
    atom_ty3     text,
    atom_ty4     text,
    force_const  float,
    n            float,
    delta        float
  );

  execute 'copy ImproperParams ' ||
          ' from ' || quote_literal(file_prefix || '_imp_params.csv') ||
          ' with csv header';

  drop table if exists NonBondedParams;
  create temp table NonBondedParams (
      atom_type  text,
      rmin       float,
      eps        float
  );

  execute 'copy NonBondedParams ' ||
          ' from ' || quote_literal(file_prefix || '_nb_params.csv') ||
          ' with csv header';



end
$$ language plpgsql;

drop function if exists load_config_atoms(bigint);
create or replace function load_config_atoms(c_id bigint)
returns void as $$
begin
  insert into AtomTypeIDs (atom_type)
  select distinct atom_type
  from AtomParams
  where not exists (select '0' 
                    from AtomTypeIDs 
                    where AtomTypeIDs.atom_type = AtomParams.atom_type);

  insert into ConfigAtoms (config_id, residue_name, atom_type_id, atom_name, charge)
  select c_id as config_id, residue_name, atom_type_id, atom_name, charge
  from AtomParams A, AtomTypeIDs T
  where A.atom_type = T.atom_type;


end
$$ language plpgsql;

drop function if exists load_config_bonds(bigint);
create or replace function load_config_bonds(c_id bigint)
returns void as $$
begin
  insert into ConfigBonds
  select c_id as config_id, 
         A1.atom_type_id as atom_type_id1, A2.atom_type_id as atom_type_id2, 
         bond_const, bond_length
  from BondParams B, AtomTypeIDs A1, AtomTypeIDs A2
  where    (B.atom_ty1 = A1.atom_type and B.atom_ty2 = A2.atom_type)
        or (B.atom_ty2 = A1.atom_type and B.atom_ty1 = A2.atom_type);

end
$$ language plpgsql;

drop function if exists load_config_angles(bigint);
create or replace function load_config_angles(c_id bigint)
returns void as $$
begin
  insert into ConfigAngles
  select c_id as config_id,
         A1.atom_type_id as atom_type_id1, 
         A2.atom_type_id as atom_type_id2, 
         A3.atom_type_id as atom_type_id3,
         NP.angle_const, NP.angle
  from AtomTypeIDs A1, AtomTypeIDs A2, AtomTypeIDs A3,
       AngleParams NP
  where      (A1.atom_type = NP.atom_ty1 and A2.atom_type = NP.atom_ty2
                                         and A3.atom_type = NP.atom_ty3)
          or (A3.atom_type = NP.atom_ty1 and A2.atom_type = NP.atom_ty2
                                         and A1.atom_type = NP.atom_ty3);

end
$$ language plpgsql;

drop function if exists load_config_dihedrals(bigint);
create or replace function load_config_dihedrals(c_id bigint)
returns void as $$
declare
  rec record;
begin



  drop table if exists ConfigDihedralsTemp;
  create temp table ConfigDihedralsTemp (
    atom_type_id1  int,
    atom_type_id2  int,
    atom_type_id3  int,
    atom_type_id4  int,
    improper    boolean,
    force_const float,
    n           float,
    delta       float,
    primary key  (atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4)
  );



  insert into ConfigDihedralsTemp
  select distinct on (atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4) * from
   ((select distinct on (atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4)
           A1.atom_type_id as atom_type_id1, A2.atom_type_id as atom_type_id2,
           A3.atom_type_id as atom_type_id3, A4.atom_type_id as atom_type_id4,
           false as improper,
           DP.force_const, DP.n, DP.delta
     from AtomTypeIDs A1, AtomTypeIDs A2, AtomTypeIDs A3, AtomTypeIDs A4,
          DihedralParams DP
        where
           ( (A1.atom_type = DP.atom_ty1 or  DP.atom_ty1 = 'X')
         and (A4.atom_type = DP.atom_ty4 or  DP.atom_ty4 = 'X')
         and (A2.atom_type = DP.atom_ty2 or  DP.atom_ty2=  'X')
         and  A3.atom_type = DP.atom_ty3 )
    )
    union
    (select distinct on (atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4)
           A1.atom_type_id as atom_type_id1, A2.atom_type_id as atom_type_id2,
           A3.atom_type_id as atom_type_id3, A4.atom_type_id as atom_type_id4,
           false as improper,
           DP.force_const, DP.n, DP.delta
     from AtomTypeIDs A1, AtomTypeIDs A2, AtomTypeIDs A3, AtomTypeIDs A4,
          DihedralParams DP
        where
           ( (A1.atom_type = DP.atom_ty1 or  DP.atom_ty1 = 'X')
         and (A4.atom_type = DP.atom_ty4 or  DP.atom_ty4 = 'X')
         and (A3.atom_type = DP.atom_ty2 or  DP.atom_ty2 = 'X')
         and  A2.atom_type = DP.atom_ty3 )
    )
    union
    (select distinct on (atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4)
           A1.atom_type_id as atom_type_id1, A2.atom_type_id as atom_type_id2,
           A3.atom_type_id as atom_type_id3, A4.atom_type_id as atom_type_id4,
           false as improper,
           DP.force_const, DP.n, DP.delta
     from AtomTypeIDs A1, AtomTypeIDs A2, AtomTypeIDs A3, AtomTypeIDs A4,
          DihedralParams DP
        where
           ( (A4.atom_type = DP.atom_ty1 or  DP.atom_ty1 = 'X')
         and (A1.atom_type = DP.atom_ty4 or  DP.atom_ty4 = 'X')
         and (A2.atom_type = DP.atom_ty2 or  DP.atom_ty2=  'X')
         and  A3.atom_type = DP.atom_ty3 )
    )
    union
    (select distinct on (atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4)
           A1.atom_type_id as atom_type_id1, A2.atom_type_id as atom_type_id2,
           A3.atom_type_id as atom_type_id3, A4.atom_type_id as atom_type_id4,
           false as improper,
           DP.force_const, DP.n, DP.delta
     from AtomTypeIDs A1, AtomTypeIDs A2, AtomTypeIDs A3, AtomTypeIDs A4,
          DihedralParams DP
        where
           ( (A4.atom_type = DP.atom_ty1 or  DP.atom_ty1 = 'X')
         and (A1.atom_type = DP.atom_ty4 or  DP.atom_ty4 = 'X')
         and (A3.atom_type = DP.atom_ty2 or  DP.atom_ty2=  'X')
         and  A2.atom_type = DP.atom_ty3 )
    ))a;




  drop table if exists ConfigImproperTemp;
  create temp table ConfigImproperTemp (
    atom_type_id1  int,
    atom_type_id2  int,
    atom_type_id3  int,
    atom_type_id4  int,
    improper    boolean,
    force_const float,
    n           float,
    delta       float,
    primary key  (atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4)
  );

  insert into ConfigImproperTemp
  select distinct on (atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4) * from
    ((select distinct on (atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4)
           A1.atom_type_id as atom_type_id1, A2.atom_type_id as atom_type_id2,
           A3.atom_type_id as atom_type_id3, A4.atom_type_id as atom_type_id4,
           true as improper,
           DP.force_const, DP.n, DP.delta
     from AtomTypeIDs A1, AtomTypeIDs A2, AtomTypeIDs A3, AtomTypeIDs A4,
          ImproperParams DP
        where
           ( (A1.atom_type = DP.atom_ty1 or  DP.atom_ty1 = 'X')
         and (A4.atom_type = DP.atom_ty4 or  DP.atom_ty4 = 'X')
         and (A2.atom_type = DP.atom_ty2 or  DP.atom_ty2=  'X')
         and  A3.atom_type = DP.atom_ty3 )
    )
    union
    (select distinct on (atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4)
           A1.atom_type_id as atom_type_id1, A2.atom_type_id as atom_type_id2,
           A3.atom_type_id as atom_type_id3, A4.atom_type_id as atom_type_id4,
           true as improper,
           DP.force_const, DP.n, DP.delta
     from AtomTypeIDs A1, AtomTypeIDs A2, AtomTypeIDs A3, AtomTypeIDs A4,
          ImproperParams DP
        where
           ( (A1.atom_type = DP.atom_ty1 or  DP.atom_ty1 = 'X')
         and (A4.atom_type = DP.atom_ty4 or  DP.atom_ty4 = 'X')
         and (A3.atom_type = DP.atom_ty2 or  DP.atom_ty2 = 'X')
         and  A2.atom_type = DP.atom_ty3 )
    )
    union
    (select distinct on (atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4)
           A1.atom_type_id as atom_type_id1, A2.atom_type_id as atom_type_id2,
           A3.atom_type_id as atom_type_id3, A4.atom_type_id as atom_type_id4,
           true as improper,
           DP.force_const, DP.n, DP.delta
     from AtomTypeIDs A1, AtomTypeIDs A2, AtomTypeIDs A3, AtomTypeIDs A4,
          ImproperParams DP
        where
           ( (A4.atom_type = DP.atom_ty1 or  DP.atom_ty1 = 'X')
         and (A1.atom_type = DP.atom_ty4 or  DP.atom_ty4 = 'X')
         and (A2.atom_type = DP.atom_ty2 or  DP.atom_ty2=  'X')
         and  A3.atom_type = DP.atom_ty3 )
    )
    union
    (select distinct on (atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4)
           A1.atom_type_id as atom_type_id1, A2.atom_type_id as atom_type_id2,
           A3.atom_type_id as atom_type_id3, A4.atom_type_id as atom_type_id4,
           true as improper,
           DP.force_const, DP.n, DP.delta
     from AtomTypeIDs A1, AtomTypeIDs A2, AtomTypeIDs A3, AtomTypeIDs A4,
          ImproperParams DP
        where
           ( (A4.atom_type = DP.atom_ty1 or  DP.atom_ty1 = 'X')
         and (A1.atom_type = DP.atom_ty4 or  DP.atom_ty4 = 'X')
         and (A3.atom_type = DP.atom_ty2 or  DP.atom_ty2=  'X')
         and  A2.atom_type = DP.atom_ty3 )
    ))a;

--   (select 
--            A1.atom_type_id as atom_type_id1, A2.atom_type_id as atom_type_id2,
--            A3.atom_type_id as atom_type_id3, A4.atom_type_id as atom_type_id4,
--            true as improper,
--            DP.force_const, DP.delta
--     from AtomTypeIDs A1, AtomTypeIDs A2, AtomTypeIDs A3, AtomTypeIDs A4,
--          ImproperParams DP
--     where
--           ( (A2.atom_type = DP.atom_ty2 or  DP.atom_ty2 = 'X')
--         and (A3.atom_type = DP.atom_ty3 or  DP.atom_ty3 = 'X')
--         and  A1.atom_type = DP.atom_ty1 and A4.atom_type = DP.atom_ty4)
--       or
--           ( (A2.atom_type = DP.atom_ty2 or  DP.atom_ty2 = 'X')
--         and (A3.atom_type = DP.atom_ty3 or  DP.atom_ty3 = 'X')
--         and  A4.atom_type = DP.atom_ty1 and A1.atom_type = DP.atom_ty4)
--       or
--           ( (A3.atom_type = DP.atom_ty2 or  DP.atom_ty2 = 'X')
--         and (A2.atom_type = DP.atom_ty3 or  DP.atom_ty3 = 'X')
--         and  A1.atom_type = DP.atom_ty1 and A4.atom_type = DP.atom_ty4 )
--       or
--           ( (A3.atom_type = DP.atom_ty2 or  DP.atom_ty2 = 'X')
--         and (A2.atom_type = DP.atom_ty2 or  DP.atom_ty3 = 'X')
--         and  A4.atom_type = DP.atom_ty1 and A1.atom_type = DP.atom_ty4 )
--     )a;
 
  insert into ConfigDihedrals
  select distinct on (atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4) c_id,* 
  from ((select * from ConfigDihedralsTemp) union 
        (select * from ConfigImproperTemp))a 
  order by atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4,improper desc;

end
$$ language plpgsql;

drop function if exists load_config_nonbondedpairs(bigint);
create or replace function load_config_nonbondedpairs(c_id bigint)
returns void as $$
begin
  create temp table ConfigNonBondedTemp (
    atom_type_id  int,
    atom_type     text,
    eps           float,
    rmin          float
  );

  insert into ConfigNonBondedTemp
  select distinct on (atom_type) atom_type_id, atom_type, eps, rmin from
    (select A.atom_type_id, A.atom_type, NBP.eps, NBP.rmin, 
            A.atom_type = NBP.atom_type_pattern as exact
     from (
        select replace(replace(replace(replace(
                atom_type,
                '%', '_'),
                '*', '%'),
                '#', '[0-9]*'),
                '+', '[0-9]')
              as atom_type_pattern, eps, rmin
       from NonBondedParams
     ) NBP, AtomTypeIDs A
     where 
         (A.atom_type    = NBP.atom_type_pattern or
          A.atom_type like NBP.atom_type_pattern))a
  order by atom_type, exact desc;


  insert into ConfigVDWNonBonded
  select c_id                                                             as config_id,
         NBT.atom_type_id                                                 as atom_type_id1, 
         NBT1.atom_type_id                                                as atom_type_id2,
         (NBT.rmin + NBT1.rmin)                                           as vdw_rmin,  
         SQRT(NBT.eps * NBT1.eps)                                         as vdw_eps
  from ConfigNonBondedTemp NBT, ConfigNonBondedTemp NBT1
  where NBT.rmin <> 0;


end
$$ language plpgsql;



----------------------------------------------------------------------------------------------
-- Protein Loader
----------------------------------------------------------------------------------------------

drop function if exists load_protein(bigint, text, text);

create or replace function load_protein(
  c_id        bigint,
  p_id        bigint,
  file_prefix text
) returns void as $$
declare
  val          int;
  rowrec       record;
begin


  --------------------------------------------
  -- Atoms for a specific protein-config pair

  create temp table ProteinAtomsTemp (
    atom_id      int,
    atom_type    text,
    atom_name    text,
    residue_id   int,
    residue_name text,
    segment_name text
  );

  execute 'copy ProteinAtomsTemp ' ||
          ' from ' || quote_literal(file_prefix || '_atom_ids.csv') ||
          ' with csv header';

  insert into AtomTypeIDs (atom_type)
  select distinct atom_type
  from ProteinAtomsTemp
  where not exists (select '0'
                    from AtomTypeIDs
                    where AtomTypeIDs.atom_type = ProteinAtomsTemp.atom_type);



  insert into ProteinAtoms (config_id, protein_id, atom_id, atom_type,
                            atom_type_id,  atom_name, residue_id, 
                            residue_name, segment_name)
  select c_id, p_id, atom_id, T.atom_type,
         atom_type_id, atom_name, residue_id, 
         residue_name, segment_name 
  from ProteinAtomsTemp T, AtomTypeIDs A
  where T.atom_type = A.atom_type;

  insert into ProteinSizes
  select c_id, p_id, count(*) from ProteinAtomsTemp;

  drop table ProteinAtomsTemp;


  --------------------------------------------
  -- Bonds for a specific protein-config pair

  create temp table ProteinBondIds (
    atom_id1     int,
    atom_id2     int
  );


  execute 'copy ProteinBondIds ' ||
          ' from ' || quote_literal(file_prefix || '_bond_ids.csv') ||
          ' with csv header';

  insert into ProteinBonds
  select c_id, p_id,
         A1.atom_id as atom_id1, A2.atom_id as atom_id2
  from  ProteinAtoms A1, ProteinAtoms A2, ProteinBondIds B
  where
         A1.config_id = c_id and A1.protein_id = p_id and 
         A2.config_id = c_id and A2.protein_id = p_id and 
         A1.atom_id = B.atom_id1 and A2.atom_id = B.atom_id2;

  drop table ProteinBondIds;

  ---------------------------------------------
  -- Angles for a specific protein-config pair

  create temp table ProteinAngleIds (
    atom_id1     int,
    atom_id2     int,
    atom_id3     int
  );

  execute 'copy ProteinAngleIds ' ||
          ' from ' || quote_literal(file_prefix || '_angle_ids.csv') ||
          ' with csv header';

  insert into ProteinAngles
  select c_id, p_id,
         A1.atom_id as atom_id1, A2.atom_id as atom_id2, A3.atom_id as atom_id3
  from   ProteinAtoms A1, ProteinAtoms A2, ProteinAtoms A3, ProteinAngleIds N
  where
         A1.config_id = c_id and 
         A1.protein_id = p_id and 
         A2.config_id = c_id and 
         A2.protein_id = p_id and 
         A3.config_id = c_id and 
         A3.protein_id = p_id and 
         A1.atom_id = N.atom_id1 and A2.atom_id = N.atom_id2 and A3.atom_id = N.atom_id3;

  drop table ProteinAngleIds;


  ------------------------------------------------
  -- Dihedrals for a specific protein-config pair

  create temp table ProteinDihedralIds (
    atom_id1     int,
    atom_id2     int,
    atom_id3     int,
    atom_id4     int
  );

  execute 'copy ProteinDihedralIds ' ||
          ' from ' || quote_literal(file_prefix || '_di_ids.csv') ||
          ' with csv header';

  insert into ProteinDihedrals
  select c_id, p_id, atom_id1, atom_id2, atom_id3, atom_id4
  from ProteinDihedralIds;


--  select c_id, p_id, 
--         A1.atom_id as atom_id1, A2.atom_id as atom_id2,
--         A3.atom_id as atom_id3, A4.atom_id as atom_id4,
--         CD.improper
--  from ProteinAtoms A1, ProteinAtoms A2, ProteinAtoms A3, ProteinAtoms A4,
--       ProteinDihedralIds D, ConfigDihedrals CD
--  where  
--         A1.config_id = c_id and A1.protein_id = p_id and 
--         A2.config_id = c_id and A2.protein_id = p_id and 
--         A3.config_id = c_id and A3.protein_id = p_id and 
--         A4.config_id = c_id and A4.protein_id = p_id and 
--         A1.atom_id = D.atom_id1 and A2.atom_id = D.atom_id2 and
--         A3.atom_id = D.atom_id3 and A4.atom_id = D.atom_id4 and
--         CD.config_id = c_id and
--         A1.atom_type_id = CD.atom_type_id1 and A2.atom_type_id = CD.atom_type_id2 and
--         A3.atom_type_id = CD.atom_type_id3 and A4.atom_type_id = CD.atom_type_id4;

  drop table ProteinDihedralIds;

  drop table if exists TempFeatures;
  create temp table TempFeatures (
    config_id         bigint,
    protein_id        bigint,
    atom_id1          int,
    atom_id2          int,
    atom_id3          int,
    atom_id4          int,
    feature_id        int
  );
  
  insert into TempFeatures (config_id, protein_id, atom_id1, atom_id2, atom_id3, atom_id4, feature_id) 
  select c_id, p_id, D.atom_id1, D.atom_id2, D.atom_id3, D.atom_id4, row_number() over (order by D.atom_id1, D.atom_id2, D.atom_id3, D.atom_id4)
  from ProteinDihedrals D,
       ProteinAtoms M1, ProteinAtoms M2, ProteinAtoms M3, ProteinAtoms M4
  where (D.atom_id1 = M1.atom_id)
  and   (D.atom_id2 = M2.atom_id)
  and   (D.atom_id3 = M3.atom_id)
  and   (D.atom_id4 = M4.atom_id)
  and   D.protein_id = p_id and   D.config_id = c_id
  and   M1.protein_id = p_id and   M1.config_id = c_id
  and   M2.protein_id = p_id and   M2.config_id = c_id
  and   M3.protein_id = p_id and   M3.config_id = c_id
  and   M4.protein_id = p_id and   M4.config_id = c_id
  and ((M1.atom_name = 'N' and M2.atom_name = 'CA' and M3.atom_name = 'C' and M4.atom_name = 'N')
    or (M1.atom_name = 'C' and M2.atom_name = 'N' and M3.atom_name = 'CA' and M4.atom_name = 'C'))
  -- Ordering the dihedrals according to the residue_id's conforming with the 
  -- simulator input templates
  order by M1.residue_id, M2.residue_id, M3.residue_id, M4.residue_id;


  insert into ProteinFeatures (config_id, protein_id, atom_id1, atom_id2, atom_id3, atom_id4, feature_id)
  select config_id, protein_id, atom_id1, atom_id2, atom_id3, atom_id4, feature_id from TempFeatures;

end;
$$ language plpgsql;

insert into AminoAcids values ('A','ALA','H');
insert into AminoAcids values ('R','ARG','C');
insert into AminoAcids values ('N','ASN','P');
insert into AminoAcids values ('D','ASP','C');
insert into AminoAcids values ('C','CYS','P');
insert into AminoAcids values ('E','GLU','C');
insert into AminoAcids values ('Q','GLN','P');
insert into AminoAcids values ('G','GLY','H');
insert into AminoAcids values ('H','HSD','P');
insert into AminoAcids values ('I','ILE','H');
insert into AminoAcids values ('L','LEU','H');
insert into AminoAcids values ('K','LYS','C');
insert into AminoAcids values ('M','MET','P');
insert into AminoAcids values ('F','PHE','P');
insert into AminoAcids values ('P','PRO','P');
insert into AminoAcids values ('S','SER','P');
insert into AminoAcids values ('T','THR','P');
insert into AminoAcids values ('W','TRP','P');
insert into AminoAcids values ('Y','TYR','P');
insert into AminoAcids values ('V','VAL','H');


insert into ClusterConfig values ('localhost',1,'straight_md');
--insert into ClusterConfig values ('qp6',10,'straight_md');
--insert into ClusterConfig values ('mddb',20,'md_inv');


drop function if exists generate_actions(n int);

create or replace function generate_actions(n int) returns void as $$
declare
  r float;
  i int;
begin
  for i in (select * from generate_series(0, n, 1)) loop
    r = i/cast(n as float);
    insert into RL_Actions (peak_ratio) values (r);
  end loop;
end;
$$ language plpgsql;

select generate_actions(5);

insert into ControlParams (param_name, param_value) 
values ('histo_res_hi', '90'),
       ('histo_res_lo', '30'),
       ('num_partitions', '2'),
       ('cache_size', '6');    -- Cache Size in total





