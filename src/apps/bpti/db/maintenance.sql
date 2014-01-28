-----------------------------------------------------------------------------------
-- Generic Trajectory loader
-----------------------------------------------------------------------------------


-- Worker-side Trajectory Loader
drop function if exists create_NewAtomPositions_table(text, int);
create or replace function create_NewAtomPositions_table(
  t_name text,
  t_id int
)
returns void as $$
begin
  execute 'drop table if exists ' || t_name || ';';
  execute 'create temp table ' || t_name || '(
             trj_id      bigint,
             t           int,
             atom_id     int,
             x           float,
             y           float,
             z           float,
             primary key (trj_id, t, atom_id)
          );';
  --update Trajectories set (status, table_name) = ('Loaded', t_name) where trj_id = t_id;
end;
$$ language plpgsql;

drop function if exists create_NewEnergies_table(text);
create or replace function create_NewEnergies_table(
  t_name text
)
returns void as $$
begin
  execute 'drop table if exists ' || t_name || ';';
  execute 'create temp table ' || t_name || ' (
             trj_id        bigint,
             t             int,
             bond          float,
             angle         float,
             dihed         float,
             imprp         float,
             elect         float,
             vdw           float,
             boundary      float,
             misc          float,
             kinetic       float,
             total         float,
             temp          float,
             potential     float,
             total3        float,
             tempavg       float
          );';
  --update Trajectories set (status, table_name) = ('Loaded', t_name) where trj_id = t_id;
end;
$$ language plpgsql;




drop function if exists trajectories_insert(e_id bigint, bf int, iter int, start_from_trj int, tran_id int);
create or replace function trajectories_insert(e_id bigint, bf int, iter int, start_from_trj int, transition_id int)
returns integer as $$
declare
  curr_trj_id integer;
begin
  lock table Trajectories in exclusive mode;

  insert into Trajectories (expjob_id, biased, biased_feature, iter_id, prev_trj_id, tran_id) 
  values (e_id, bf is not null, bf, iter, start_from_trj, transition_id);

  select max(trj_id) into curr_trj_id from Trajectories;

  return curr_trj_id;
end;
$$ language plpgsql;




drop function if exists compute_features(bigint, bigint, text, text);
create or replace function compute_features(
  e_id bigint,
  ap_t_name text,
  pf_t_name text
) returns void as $$
declare
  c_id bigint;
  p_id bigint;
  st text;
begin

  select config_id, protein_id into c_id, p_id
  from Subspaces S, ExpJobs E
  where E.expjob_id = e_id and
        S.subspace_id = E.subspace_id;


  execute 'drop table if exists ' || pf_t_name;
  execute 'create temp table ' || pf_t_name || ' as 
           select P1.trj_id, P1.t,
                  BD.feature_id,
                  dihedral_angle(P1.x,P1.y,P1.z,
                                 P2.x,P2.y,P2.z,
                                 P3.x,P3.y,P3.z,
                                 P4.x,P4.y,P4.z) as phi_psi
           from ProteinFeatures BD, ' || '
               ' || ap_t_name || ' P1, ' || ap_t_name || ' P2,
               ' || ap_t_name || ' P3, ' || ap_t_name || ' P4
           where P1.t = P2.t           and P1.t = P3.t           and P1.t = P4.t
           and   P1.trj_id = P2.trj_id and P1.trj_id = P3.trj_id and P1.trj_id = P4.trj_id
           and   BD.config_id = ' || c_id || 'and BD.protein_id = ' || p_id || '
           and   BD.atom_id1 = P1.atom_id 
           and   BD.atom_id2 = P2.atom_id
           and   BD.atom_id3 = P3.atom_id 
           and   BD.atom_id4 = P4.atom_id';

  --execute 'drop table if exists PhiPsiTemp';
  --execute 'create temp table PhiPsiTemp as select * from ' || pf_t_name;
end;
$$ language plpgsql;



drop function if exists create_NewHistogram(pf text, hs text);
create or replace function create_NewHistogram(pf text, hs text) returns void as $$
declare
  num_features int;
  ss_id int;
  d_id int;
begin


  execute 'drop table if exists ' || hs;
  execute 'create temp table ' || hs || ' as
           select trj_id,
                  id_array, 
                  region_id,
                  count(*) as n,
                  sum(case when is_first_point then 1 else 0 end) as s
           from (
                 select S.trj_id, 
                        array_agg(bucket_id order by B.feature_id) as id_array,
                        array_agg(bucket_id/(H.buckets_per_dim[B.feature_id]/H.regions_per_dim[B.feature_id]) 
                                  order by B.feature_id) as region_id,
                        bool_and(S.t = 1) as is_first_point
                 from  (select * from ' || pf || ' 
                       )S, 
                       Buckets B, HistogramSpec H
                 where H.spec_id = 1
                 and   B.num_dims = H.num_dims
                 and   S.feature_id = B.feature_id
                 and   S.phi_psi between B.bucket_start and B.bucket_end
                 group by S.trj_id, S.t
           ) as R
           group by R.trj_id, R.id_array, R.region_id';

end
$$ language plpgsql;



drop function if exists process_new_data();
create or replace function process_new_data(
  e_id bigint
) returns void as $$
declare
  st text;
  sub_queries text[];
  hs text;
begin
  --select array_agg('(select * from ' || ap_t_name || ')')  into sub_queries
  --from tables_to_load;
  --st = E'(' || array_to_string(sub_queries, E'\n union all\n') || E')a;';

  --execute 'insert into AtomPositions select * from ' || st;

  --select array_agg('(select * from ' || pf_t_name || ')')  into sub_queries
  --from tables_to_load;
  --st = E'(' || array_to_string(sub_queries, E'\n union all\n') || E')a;';

  --execute 'insert into ConformationSpace select * from ' || st;

  select array_agg('(select * from ' || hs_t_name || ')')  into sub_queries
  from tables_to_load;

  st = E'(' || array_to_string(sub_queries, E'\n union all\n') || E')';

  --execute 'select subspace_id, desc_id, id_array, region_id, 
  --                max(point_id_pair) as point_id_pair,
  --                sum(n) as n,
  --                sum(s) as s
  --         from (' || st || ')a
  --         group by subspace_id, desc_id, id_array, region_id';

  perform update_Histogram(e_id, st);
end;
$$ language plpgsql;


drop function if exists update_Histogram(e_id bigint, st text);
create or replace function update_Histogram(
  e_id bigint,
  st text
) returns void as $$
declare
  rec record;
  stats record;
  total_count int;
  c int = 0;
  ss_id int;
  d_id int;
begin

  --lock table HistogramStats in exclusive mode;
  --lock table Histogram in exclusive mode;
  --lock table HistogramRegionScores in exclusive mode;

  select subspace_id, desc_id into ss_id, d_id
  from ExpJobs E
  where E.expjob_id = e_id; 

  for rec in execute (st) loop
    update Histogram set n = rec.n + n, s = rec.s + s 
      where subspace_id = ss_id and
            desc_id = d_id and
            id_array = rec.id_array;

    if not found then
      insert into Histogram (subspace_id, desc_id, id_array, region_id, n, s, trj_id, t)
        values (ss_id, d_id, rec.id_array, rec.region_id, 
                rec.n, rec.s, rec.point_id_pair[1], rec.point_id_pair[2]);
      c = c+1;
    end if;

    update HistogramRegionScores 
    set n     = n + rec.n, 
        s     = s + rec.s, 
        score = (n + rec.n) * (s + rec.s + 1) 
    where subspace_id = ss_id and
          desc_id = d_id and
          region_id  = rec.region_id;


    if not found then
      insert into HistogramRegionScores (subspace_id, desc_id, region_id, n, s, score) 
        values (ss_id, d_id, rec.region_id, rec.n, rec.s, rec.n* (rec.s+1));
    end if;
  end loop;

--  select * into stats from HistogramStats 
--  where subspace_id = ss_id and desc_id = d_id  order by id desc limit 1;
--  if not found then
--    insert into HistogramStats (subspace_id, desc_id, curr_size, max_order_of_mag)
--    select ss_id, d_id, 0, sum(ln(num_buckets))
--    from (select unnest(buckets_per_dim) as num_buckets
--         from HistogramSpec where spec_id = 1)a;
--    select * into stats from HistogramStats order by id desc limit 1;
--  else
--    insert into HistogramStats (subspace_id, desc_id, curr_size, max_order_of_mag) 
--    values (ss_id, d_id, stats.curr_size + c, stats.max_order_of_mag);
--  end if;

end
$$ language plpgsql;

create or replace function copy_histograms(ss_id int, d_id int, x int, y int) returns void as $$
declare
begin

  drop table if exists temp_histo;
  create temp table temp_histo as
  select id_array[x] as a, id_array[y] as b, sum(n) as n, sum(s) as s
  from Histogram  H, Trajectories T, ExpJobs E
  where H.trj_id = T.trj_id and
        T.expjob_id = E.expjob_id and
        E.subspace_id = ss_id and 
        E.desc_id = d_id
  group by id_array[x], id_array[y];

  copy temp_histo to '/tmp/histogram.csv' with delimiter ',';

  drop table if exists temp_histo;
  create temp table temp_histo as
  select region_id[x] as a, region_id[y] as b, sum(n) as n, sum(s) as s, sum(score) as score
  from HistogramRegionScores  H, Trajectories T, ExpJobs E
  where H.trj_id = T.trj_id and
        T.expjob_id = E.expjob_id and
        E.subspace_id = ss_id and 
        E.desc_id = d_id
  group by region_id[x], region_id[y];

  copy temp_histo to '/tmp/histo_scores.csv' with delimiter ',';
end;
$$ language plpgsql;





create or replace function copy_histograms(ss_id int, d_id int, x int, y int, z int) returns void as $$
declare
begin

  drop table if exists temp_histo;
  create temp table temp_histo as
  select id_array[x] as a, id_array[y] as b, id_array[z] as c, sum(n) as n, sum(s) as s
  from Histograms H, Trajectories T, ExpJobs E
  where H.trj_id = T.trj_id and
        T.expjob_id = E.expjob_id and
        E.subspace_id = ss_id and 
        E.desc_id = d_id
  group by id_array[x], id_array[y], id_array[z];

  copy temp_histo to '/tmp/histogram3d.csv' with delimiter ',';

  drop table if exists temp_histo;
  create temp table temp_histo as
  select region_id[x] as a, region_id[y] as b, region_id[z] as c, sum(n) as n, sum(s) as s, sum(score) as score
  from HistogramRegionScores H, Trajectories T, ExpJobs E
  where H.trj_id = T.trj_id and
        T.expjob_id = E.expjob_id and
        E.subspace_id = ss_id and 
        E.desc_id = d_id
  group by region_id[x], region_id[y], region_id[z];

  copy temp_histo to '/tmp/histo_scores3d.csv' with delimiter ',';
end;
$$ language plpgsql;




create or replace function copy_histograms() returns void as $$
declare
begin

  copy Histogram to '/tmp/histogram.csv' with delimiter ',';

  copy HistogramRegionScores to '/tmp/histo_scores.csv' with delimiter ',';
end;
$$ language plpgsql;



