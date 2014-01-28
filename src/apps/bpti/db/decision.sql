-- watch -n 2 "psql -X -q  -1 -v ON_ERROR_STOP=0 --pset pager=off -d test_db -c 'select count(*) from FinishedTrajectories'"

--  watch -n 2 "psql -X -q  -1 -v ON_ERROR_STOP=0 --pset pager=off -d test_db -c 'select count(*) from gmqueue'"
-- 



create or replace function compute_bumpiness(iter bigint) returns float as $$
declare
  res               float;
  bins_per_region   int;
  ss_id             int;
  histo_res_hi int;
begin
  --perform cache_load(iter, 'r', 'compute_bumpiness');

  select sqrt(variance(R.n))/avg(R.n) into res
  from
    (select log(n+1) as n
     from TempIterHistograms H
     where H.iter_id = iter 
    ) R;

  return res;
end
$$ language plpgsql;


drop function if exists compute_bumpiness2(bigint);
create or replace function compute_bumpiness2(iter bigint) returns float as $$
declare
  res               float=0;
  num_bins          int = 15;
  bins_per_region   int;
  ss_id             int;
  histo_res_hi int;
  num_features      int;
  j                 int;
  b                 float;
begin
  select cast(param_value as int) into histo_res_hi
  from ControlParams where param_name = 'histo_res_hi';

  bins_per_region = histo_res_hi/num_bins;

  perform cache_load(iter, 'r', 'compute_bumpiness');

  select subspace_id into ss_id from RL_History R where R.iter_id = iter;


  select count(*) into num_features
  from ProteinFeatures P, Subspaces S
  where S.subspace_id = ss_id and
        P.config_id = S.config_id and
        P.protein_id = S.protein_id;

  for j in (select generate_series(1,num_features-1))
  loop
    select sqrt(variance(R.n))/avg(R.n) into b
    from
      (select bucket_id[j] as key1, bucket_id[j+1] as key2,
              log(sum(n)) as n
       from IterHistogramsCache H
       where H.iter_id = iter
       group by key1,key2
      ) R;
  
    res = res+b;
  end loop;

  return res/(num_features-1);
end
$$ language plpgsql;



create or replace function compare_histograms(curr_iter bigint)
returns boolean
as $$
declare
  res boolean;
begin
  create temp table HistoA as
  select id_array as key, sum(n) as n
  from Trajectories T, Histograms B
  where T.trj_id = B.trj_id and
        T.iter_id <= curr_iter
  group by id_array;

  create temp table HistoB as
  select array(select a/2 from unnest(bucket_id)a) as key, 
         sum(n) as n
  from IterHistograms
  where iter_id = curr_iter
  group by key;


  select bool_and(A.n - B.n = 0) into res
  from HistoA A full outer join HistoB B
  on A.key = B.key;

  return res;
end
$$ language plpgsql;



create or replace function update_histograms(curr_iter bigint)
returns void
as $$
declare
  histo_res_hi int;
  histo_res_lo int;
  num_partitions int;
  ss_id int;
  prev_iter int;
  last_conf_id int;
  nf int;
  np int;
  nb int;
  h2l_ratio int;
  max_count int;
  min_count int;
  threshold int;
  num_dims  int;
  histo_size int;
  histo_size_limit int = 20000;
begin
  select subspace_id into ss_id from RL_History R where R.iter_id = curr_iter;

  select count(*) into num_dims
  from ProteinFeatures P, Subspaces S 
  where P.config_id   = S.config_id and 
        P.protein_id  = S.protein_id and
        S.subspace_id = ss_id;

  select iter_id into prev_iter 
  from RL_History R 
  where R.iter_id < curr_iter and
        R.subspace_id = ss_id 
  order by R.iter_id desc limit 1;

  select cast(param_value as int) into histo_res_hi
  from ControlParams where param_name = 'histo_res_hi';

  select cast(param_value as int) into histo_res_lo
  from ControlParams where param_name = 'histo_res_lo';

  h2l_ratio = histo_res_hi/histo_res_lo;

  select cast(param_value as int) into num_partitions
  from ControlParams where param_name = 'num_partitions';

  select conf_id into last_conf_id from HistoBinAssignments order by conf_id desc limit 1;
  if last_conf_id is null then
    last_conf_id = 0;
  end if;

  create temp table TempHistoBinAssignments as
  select row_number() over (order by C.trj_id, C.t) + last_conf_id as conf_id,
         curr_iter as iter_id,
         C.trj_id,
         C.t,
         array_agg( cast(floor(((C.phi_psi+pi())/(pi()*2)) * histo_res_hi) as int)
                       order by feature_id) as bucket_id,
         min(EN.total) as e
  from ConformationSpace C, Trajectories T, Energies EN
  where C.trj_id  = T.trj_id and
        T.iter_id = curr_iter and
        EN.trj_id  = C.trj_id and
        EN.t       = C.t
  group by C.trj_id,C.t;

  select count(*) into nf from HistoBinAssignments;

  select protein_size*nf into np from ProteinSizes P, Subspaces S
  where S.subspace_id = ss_id and
        P.config_id = S.config_id and
        P.protein_id = S.protein_id;

  insert into HistoBinAssignments (conf_id, iter_id, trj_id, t, bucket_id, e)
  select TH.conf_id, TH.iter_id, TH.trj_id, TH.t, TH.bucket_id, TH.e
  from TempHistoBinAssignments TH;

  create temp table TempNewPartitionedHistogramsHI as
  select partition_id, 
         bucket_id, 
         count(*) as n, 
         sum(case when t=1 then 1 else 0 end) as s,
         max(conf_id) as conf_id,
         sum(e) as e
  from
    (select floor(random() * num_partitions)+1 as partition_id,
            HA.*
     from TempHistoBinAssignments HA
     where HA.iter_id = curr_iter
    ) as A
  group by A.partition_id, A.bucket_id;

  create temp table TempNewPartitionedHistogram as 
  select partition_id, 
         reduce_resolution(bucket_id, h2l_ratio) as bucket_id, 
         sum(n) as n,
         sum(s) as s,
         max(conf_id) as conf_id,
         sum(e) as e
  from TempNewPartitionedHistogramsHI
  group by partition_id, bucket_id;

  perform cache_load(prev_iter, 'r', 'update_histograms (prev_iter)');

  create temp table TempIterHistograms as
  select curr_iter as iter_id,
         bucket_id,
         sum(n) as n,
         sum(s) as s,
         max(conf_id) as conf_id,
         sum(e) as e
  from (
         select bucket_id, n, s, conf_id,e from TempNewPartitionedHistogram 
         union all
         (select bucket_id, n, s, conf_id,e from IterHistogramsCache where iter_id = prev_iter
          order by e/n limit histo_size_limit)
       ) A
  group by bucket_id;
  create index TempIterHistograms_index on TempIterHistograms (bucket_id);

  create temp table TempPartitionedHistograms as
  select curr_iter as iter_id,
         A.partition_id,
         A.bucket_id,
         sum(A.n) as n,
         sum(A.s) as s,
         max(A.conf_id) as conf_id,
         sum(e) as e
  from (
         select * from TempNewPartitionedHistogram
         union all
         (select partition_id, bucket_id, n, s, conf_id, e
          from PartitionedHistogramsCache where iter_id = prev_iter
          order by e/n limit histo_size_limit*num_partitions)
       ) A
  group by A.partition_id, A.bucket_id;
  create index TempPartitionedHistograms_index on TempPartitionedHistograms (bucket_id);

  create temp table TempIterHistogramsHI as
  select curr_iter as iter_id,
         A.bucket_id,
         sum(A.n) as n,
         sum(A.s) as s,
         max(A.conf_id) as conf_id,
         sum(e) as e
  from (
         select bucket_id, n, s, conf_id, e
         from TempNewPartitionedHistogramsHI
         union all
         (select bucket_id, n, s, conf_id,e from IterHistogramsHICache where iter_id = prev_iter
          order by e/n limit histo_size_limit*h2l_ratio)
       ) A
  group by A.bucket_id;

  create index TempIterHistogramsHI_index on TempIterHistogramsHI (bucket_id);

  select count(*) into nb from TempPartitionedHistograms;
  select count(*)+nb into nb from TempIterHistograms;

  insert into SystemStats values (curr_iter, nf, np, nb);

  insert into IterHistogramsHI
  select * from TempIterHistogramsHI;

  drop table TempIterHistogramsHI;

  insert into IterHistograms
  select * from TempIterHistograms;

  insert into PartitionedHistograms
  select * from TempPartitionedHistograms;
 

  perform cache_load(curr_iter, 'w', 'update_histograms (curr_iter)');
 
end
$$ language plpgsql;

create or replace function reduce_resolution(xs int[], bins_per_group int)
returns int[]
as $$
declare
  res int[];
begin
  select array(select a/bins_per_group from unnest(xs)a) into res;

  return res;
end
$$ language plpgsql;


create or replace function compute_accuracy2(curr_iter bigint) 
returns float
as $$
declare
  res               float=0;
  min_res           float=1;
  i                 int;
  j                 int;
  cos_sim           float;
  ss_id             int;
  num_partitions    int;
  histo_res_hi int;
  bins_per_group    int;
  num_entries       int = 0;
  num_features      int;
  size_a            float;
  size_b            float;
  resolution        int = 30;
begin
  select cast(param_value as int) into histo_res_hi
  from ControlParams where param_name = 'histo_res_hi';

  bins_per_group = histo_res_hi/resolution;


  select cast(param_value as int) into num_partitions
  from ControlParams where param_name = 'num_partitions';
  perform cache_load(curr_iter, 'r', 'compute_accuracy');
  select subspace_id into ss_id from RL_History R where R.iter_id = curr_iter;
  --num_partitions = num_partitions/3;

  for i in (select generate_series(1,num_partitions)) loop

    select sqrt(sum(n*cast(n as float))) into size_a from PartitionedHistogramsCache
    where partition_id = i and iter_id = curr_iter;

    select sqrt(sum(n*cast(n as float))) into size_b 
    from (select sum(n) as n 
          from PartitionedHistogramsCache
          where partition_id != i and iter_id = curr_iter
          group by bucket_id)a;




    select sum(A.n*cast(B.n as float))/ (size_a*size_b) into cos_sim
    from (select bucket_id as key1, n
          from PartitionedHistogramsCache
          where partition_id = i and iter_id = curr_iter
         ) A,
         (select bucket_id as key1, sum(n) as n
          from PartitionedHistogramsCache
          where partition_id != i and iter_id = curr_iter
          group by bucket_id
         ) B
    where A.key1 = B.key1;


    if cos_sim is null then
      cos_sim = 0;
    end if;
    res = res + cos_sim;
  end loop;
  res = res/num_partitions;
  return res;
end
$$ language plpgsql;

create or replace function compute_accuracy(curr_iter bigint)
returns float
as $$
declare
  res               float=0;
  min_res           float=1;
  i                 int;
  j                 int;
  cos_sim           float;
  ss_id             int;
  num_partitions    int;
  histo_res_hi int;
  bins_per_group    int;
  num_entries       int = 0;
  num_features      int;
  size_a            float;
  size_b            float;
  resolution        int = 30;
begin
  select cast(param_value as int) into histo_res_hi
  from ControlParams where param_name = 'histo_res_hi';

  bins_per_group = histo_res_hi/resolution;


  select cast(param_value as int) into num_partitions
  from ControlParams where param_name = 'num_partitions';
  --perform cache_load(curr_iter, 'r', 'compute_accuracy');
  select subspace_id into ss_id from RL_History R where R.iter_id = curr_iter;
  --num_partitions = num_partitions/3;

  for i in (select generate_series(1,num_partitions)) loop

    select sqrt(sum(n*cast(n as float))) into size_a from TempPartitionedHistograms
    where partition_id = i and iter_id = curr_iter;

    select sqrt(sum(n*cast(n as float))) into size_b
    from (select sum(n) as n
          from TempPartitionedHistograms
          where partition_id != i and iter_id = curr_iter
          group by bucket_id)a;




    select sum(A.n*cast(B.n as float))/ (size_a*size_b) into cos_sim
    from (select bucket_id as key1, n
          from TempPartitionedHistograms
          where partition_id = i and iter_id = curr_iter
         ) A,
         (select bucket_id as key1, sum(n) as n
          from TempPartitionedHistograms
          where partition_id != i and iter_id = curr_iter
          group by bucket_id
         ) B
    where A.key1 = B.key1;


    if cos_sim is null then
      cos_sim = 0;
    end if;
    res = res + cos_sim;
  end loop;
  drop table TempPartitionedHistograms;
  res = res/num_partitions;
  return res;
end
$$ language plpgsql;


create or replace function compute_accuracy2(curr_iter bigint) 
returns float
as $$
declare
  res               float=0;
  min_res           float=1;
  i                 int;
  j                 int;
  cos_sim           float;
  ss_id             int;
  num_partitions    int;
  histo_res_hi int;
  bins_per_group    int;
  num_entries       int = 0;
  num_features      int;
begin
  select cast(param_value as int) into num_partitions
  from ControlParams where param_name = 'num_partitions';
  perform cache_load(curr_iter, 'r', 'compute_accuracy');
  select subspace_id into ss_id from RL_History R where R.iter_id = curr_iter;
  --num_partitions = num_partitions/3;

  select count(*) into num_features
  from ProteinFeatures P, Subspaces S
  where S.subspace_id = ss_id and
        P.config_id = S.config_id and
        P.protein_id = S.protein_id;
  --num_features = greatest(2,num_features/3);
  for j in (select generate_series(1,num_features-1))
  loop
    res = 0;
    for i in (select generate_series(1,num_partitions)) loop
      select sum(A.n*B.n)/ (sqrt(sum(A.n*A.n))* sqrt(sum(B.n*B.n))) into cos_sim
      from (select bucket_id[j] as key1, bucket_id[j+1] as key2, sum(n) as n
            from PartitionedHistogramsCache
            where partition_id = i and iter_id = curr_iter
            group by key1,key2) A 
           full outer join 
           (select bucket_id[j] as key1, bucket_id[j+1] as key2, sum(n) as n
            from IterHistogramsCache 
            where iter_id = curr_iter
            group by key1,key2) B
      on A.key1 = B.key1 and A.key2 = B.key2;
      if cos_sim is null then
        cos_sim = 0;
      end if;
      res = res + cos_sim;
    end loop;
    res = res/num_partitions;
    min_res = least(res,min_res);
  end loop;
  return min_res;
end
$$ language plpgsql;

--
--
--create or replace function compute_bumpiness(iter_id bigint) returns float as $$
--declare
--  res             float;
--  num_bins        int default 180;
--  num_groups      int default 45;
--  bins_per_group  int;
--begin
--  bins_per_group = num_bins/num_groups;
--
--  select max(c) into res
--  from
--    (select sqrt(variance(R.n))/avg(R.n) as c 
--     from
--       (select B.b, coalesce(A.n,0) as n
--        from
--          (select id_array[T.biased_feature]/bins_per_group as key, T.biased_feature, log(sum(n)) as n
--           from Trajectories T, BiasedHistograms B 
--           where T.trj_id = B.trj_id and 
--                 T.biased and
--                 T.iter_id=iter_id 
--           group by id_array[T.biased_feature]/bins_per_group, T.biased_feature
--          ) A
--        right join
--          generate_series(0,num_groups-1) B
--        on A.key = B.b
--       ) R
--     )S;
--
--  return res;
--end
--$$ language plpgsql;
--


create or replace function recompute_accuracy(ss_id int) returns void as $$
declare
  i int;
  b float;
begin
  for i in (select iter_id from RL_History where subspace_id = ss_id) loop
    select compute_accuracy(i) into b;
    update RL_History set acc_score = b where iter_id = i;
  end loop;

end
$$ language plpgsql;


create or replace function recompute_bumpiness() returns void as $$
declare
  i int;
  b float;
begin
  for i in (select iter_id from RL_History) loop
    select compute_bumpiness(i) into b;
    update RL_History set acc_score = -b where iter_id = i;
  end loop;

end
$$ language plpgsql;

drop function if exists determine_conf2sample(bigint);
create or replace function determine_conf2sample(e_id bigint) returns void as $$
declare
  ss_id       bigint;
  d_id        bigint;
  nr          int = 15;
  rs          int;
  re          int;
begin

  select subspace_id, desc_id into ss_id, d_id
  from ExpJobs E
  where E.expjob_id = e_id;

  create temp table TempHistogramRegionScores as
  select region_id, sum(n) as n, sum(s) as s, sum(n)*(sum(s)+1) as score 
  from HistogramRegionScores H, Trajectories T, ExpJobs E
  where H.trj_id = T.trj_id and
        E.expjob_id = T.expjob_id and
        E.subspace_id = ss_id and 
        E.desc_id = d_id
  group by H.region_id;

  drop table if exists Regions_to_sample;
  create temp table Regions_to_sample as
  select R.region_id, R.score + 0.1 * 
         sum(S.score/(1+tor_distl1(R.region_id, S.region_id,nr))) as score
  from (select * from TempHistogramRegionScores order by score asc limit 10) R,
       (select * from TempHistogramRegionScores order by score desc limit 10) S
  group by R.region_id, R.score;
end
$$ language plpgsql;


-- drop function if exists generate_regions_to_sample(curr_iter int, ss_id bigint);
-- create or replace function generate_regions_to_sample(curr_iter int, ss_id bigint) returns void as $$
-- declare
--   nr                int = 20;
--   min_count         int;
--   max_count         int;
--   cutoff            int;
--   ratio             float = 0.1;
--   resolution        int = 15;
--   histo_res_lo int;
--   bins_per_region   int;
--   prev_iter         int;
-- begin
--   select cast(param_value as int) into histo_res_lo
--   from ControlParams where param_name = 'histo_res_lo';
-- 
--   bins_per_region = histo_res_lo/resolution;
-- 
-- 
--   select iter_id into prev_iter
--   from RL_History R
--   where R.iter_id != curr_iter and
--         R.subspace_id = ss_id and
--         acc_score is null
--   order by R.iter_id desc limit 1;
-- 
-- 
-- 
--   --perform cache_load(prev_iter, 'r', 'generate_regions_to_sample');
-- 
--   create temp table TempHistogramRegionScores as
--   select bucket_id as region_id,
--          sum(n) as n
--   from TempIterHistograms H
--   where H.iter_id = prev_iter
--   group by region_id;
-- 
--   select max(n),min(n) into max_count, min_count from TempHistogramRegionScores;
--   cutoff = cast((ratio * (max_count - min_count) + min_count) as int);
-- 
--   drop table if exists Peaks_to_sample;
--   create temp table Peaks_to_sample as
--   select R.region_id, R.n + 0.1 *
--          sum(S.n/(1+tor_distl1(R.region_id, S.region_id,resolution))) as score
--   from (select * from TempHistogramRegionScores where n > cutoff/1000 order by n asc limit nr) R,
--        (select * from TempHistogramRegionScores where n > cutoff/1000 order by n desc limit nr) S
--   group by R.region_id, R.n;
-- 
--   drop table if exists Wells_to_sample;
--   create temp table Wells_to_sample as
--   select R.region_id, R.n + 0.1 *
--          sum(S.n/(1+tor_distl1(R.region_id, S.region_id,resolution))) as score
--   from (select * from TempHistogramRegionScores where n > cutoff order by n asc limit nr) R,
--        (select * from TempHistogramRegionScores where n > cutoff order by n desc limit nr) S
--   group by R.region_id, R.n;
-- 
-- end
-- $$ language plpgsql;
-- 

create or replace function generate_regions_to_sample(curr_iter int, ss_id bigint) returns void as $$
declare
  nr                int = 40;
  min_count         int;
  max_count         int;
  cutoff            int;
  ratio             float = 0.1;
  resolution        int;
  histo_res_lo int;
  bins_per_region   int;
begin
  --perform cache_load(prev_iter, 'r', 'generate_regions_to_sample');

  select cast(param_value as int) into resolution
  from ControlParams where param_name = 'histo_res_lo';


  select max(n),min(n) into max_count, min_count from TempIterHistograms;
  cutoff = cast((ratio * (max_count - min_count) + min_count) as int);

  drop table if exists Peaks_to_sample;
  create temp table Peaks_to_sample as
  select R.bucket_id as region_id, R.n + 0.1 *
         sum(S.n/(1+tor_distl1(R.bucket_id, S.bucket_id,resolution))) as score
  from (select * from TempIterHistograms where n > cutoff/1000 order by n asc limit nr) R,
       (select * from TempIterHistograms where n > cutoff/1000 order by n desc limit nr) S
  group by R.bucket_id, R.n;

  drop table if exists Wells_to_sample;
  create temp table Wells_to_sample as
  select R.bucket_id as region_id, R.n + 0.1 *
         sum(S.n/(1+tor_distl1(R.bucket_id, S.bucket_id,resolution))) as score
  from (select * from TempIterHistograms where n > cutoff order by n asc limit nr) R,
       (select * from TempIterHistograms where n > cutoff order by n desc limit nr) S
  group by R.bucket_id, R.n;

  drop table TempIterHistograms;
end
$$ language plpgsql;



drop table if exists RandPhiPsi;
create table RandPhiPsi (
  old_vals float[],
  new_vals float[],
  num_changes int,
  change_ratio float
);

drop table if exists StartingPoints;
create table StartingPoints (
  trj_id bigint,
  p float[]
);

drop function if exists schedule_new_jobs(bigint, text);
create or replace function schedule_new_jobs(e_id bigint, gmd_host text) 
returns int as $$
declare
   trj_id int;
   f float[];
   num_confs int;
   r_id int[];
   rpd int[];
   nd int;
   o int;
begin 
  perform determine_conf2sample(e_id);
  select count(*) from Regions_to_sample into num_confs;

  select get_random_number(1,5) into o;

  select region_id into r_id
  from Regions_to_sample 
  order by score limit 1 offset o;

  select count(*) into nd from unnest(r_id);


  select regions_per_dim into rpd 
  from HistogramSpec
  where num_dims = nd;

  select array_agg(regionid_to_rad(r_id[A.a], rpd[A.a]) order by A.a) 
  into f
  from generate_series(1,nd) A;

  update ExpJobs
  set blocked = false
  where expjob_id = e_id;


  perform queue_gearman_job(e_id, f);
  return 1;
end 
$$ language plpgsql;

drop function if exists regionid_to_rad(r_id int, num_regions int);
create or replace function regionid_to_rad(r_id int, num_regions int)
returns float as $$
begin
  return (r_id + random()) * (2*pi()/num_regions) - pi();
end 
$$ language plpgsql;

drop function if exists schedule_md_inv_jobs(ss_id int);
create or replace function schedule_md_inv_jobs(ss_id int)
returns int as $$
declare
  total_jobs int = 0;
  e_id int;
  biased_feature int;
  num_tests int = 10;
  i int;
  iter int;
begin
  select iter_id into iter
  from RL_History
  where subspace_id = ss_id and
        acc_score is null
  order by iter_id desc
  limit 1;



  select max(expjob_id) into e_id from ExpJobs E where E.subspace_id = ss_id;
  for biased_feature in (select feature_id 
                         from ProteinFeatures P, SubSpaces S 
                         where P.config_id = S.config_id and 
                               P.protein_id = S.protein_id and 
                               S.subspace_id = ss_id)
  loop
    for i in (select generate_series(1,num_tests)) loop
      perform queue_gearman_job(e_id, array[0.5], 'md_inv', biased_feature, iter);
      total_jobs = total_jobs + 1;
    end loop;
  end loop;
  return total_jobs;
end
$$ language plpgsql;


drop function if exists get_action(ss_id int, st_id int, greediness float);
create or replace function get_action(ss_id int, st_id int, greediness float)
returns int as $$
declare
  rand_val float;
  a_id int;
  num_actions int;
  po text;
begin


  select count(*) into num_actions from RL_Actions;
  a_id = num_actions;
  select policy into po from SubSpaces where subspace_id = ss_id;
  if po = 'fixed' then
    return a_id;
  end if;

  if random() < greediness then
    select Q.action_id into a_id
    from RL_QTable Q
    where Q.subspace_id = ss_id and
          Q.state_id = st_id
    order by q_value desc
    limit 1;
  else
    select action_id into a_id from rl_actions 
    order by action_id
    offset floor(random()*num_actions)
    limit 1;
  end if;

  if a_id is null then
    a_id = num_actions;
  end if;

  return a_id;
end
$$ language plpgsql;


drop function if exists find_dest_wells(ss_id int, trj_id int, initial_well int, num_wells int);
create or replace function find_dest_wells(ss_id int, trj_id int, initial_well int, num_wells int)
returns int[] as $$
declare
  dests int[];
  conf  float[];
  dest_candidate int;
  l_conf float[];
  u_conf float[];
  starting_well int;
  t int;
begin
  drop table if exists TempDestWells;
  create table TempDestWells as
  select well_id, ref_id,
         array_agg(centroid - span/2 order by feature_id) as upper,
         array_agg(centroid + span/2 order by feature_id) as lower
  from MSM_EnergyWellsPhiPsi
  where subspace_id = ss_id
  group by well_id, ref_id;
    

  starting_well = initial_well;
  for t,conf in (select C.t,array_agg(phi_psi order by feature_id)
               from ConformationSpace C 
               where C.trj_id = trj_id group by C.t order by C.t)
  loop
    for dest_candidate, l_conf, u_conf in (select well_id, upper, lower from TempDestWells) loop
      if dest_candidate != starting_well then
        if conf between l_conf and u_conf then
          starting_well = dest_candidate;
          dests = array_append(dests, dest_candidate);
          if num_wells is not NUll and array_upper(dests, 1) >= num_wells then
            return dests;
          end if;
        end if;
      end if;
    end loop;
  end loop;
  return dests;
end
$$ language plpgsql;

drop function if exists detect_transitions(ss_id int);
create or replace function detect_transitions(ss_id int)
returns void as $$
declare
  e_id              int;
  trj_id            int;
  tran_id           int;
  starting_well     int;
  dests             int[];
  start_from_trj    int;
  num_dests         int;
  d                 int;
begin
  for e_id in (select expjob_id from ExpJobs E where E.subspace_id = ss_id) loop
    select T.trj_id, M.well_id1, T.tran_id into trj_id, starting_well, tran_id
    from Trajectories T, MSM_Transitions M
    where     T.expjob_id = e_id
          and T.tran_id   = M.tran_id
    order by T.trj_id desc limit 1;

    if trj_id is not Null then
      dests = find_dest_wells(ss_id, trj_id, starting_well, Null);
      num_dests =  coalesce(array_upper(dests, 1),0);
      if num_dests > 0 then
        for d in (select unnest(dests[1:num_dests-1])) loop
          update MSM_Transitions M set well_id2 = dests[num_dests], trj_id2 = trj_id where M.tran_id = tran_id;
          insert into MSM_Transitions (subspace_id, expjob_id, well_id1, trj_id1) values (ss_id, e_id, d, trj_id);
          select M.tran_id into tran_id from MSM_Transitions M order by M.tran_id desc limit 1;
        end loop;

        update MSM_Transitions M set well_id2 = dests[num_dests], trj_id2 = trj_id where M.tran_id = tran_id;
        select M.tran_id into tran_id from MSM_Transitions M order by M.tran_id desc limit 1;
        start_from_trj = Null;
      end if;
    else
    end if;
  end loop;

end
$$ language plpgsql;




drop function if exists schedule_straight_md_jobs(ss_id int);
create or replace function schedule_straight_md_jobs(ss_id int)
returns int as $$
declare
  e_id              int;
  trj_id            int;
  tran_id           int;
  starting_well     int;
  features          float[];
  min_well_features float[];
  min_well          int;
  start_from_trj    int;
  iter              int;
  num_jobs          int = 0;
begin
  select iter_id+1 into iter from Trajectories order by iter_id desc limit 1;
  drop table if exists TempWellStats;

  create temp table TempWellStats as
  select A.well_id, coalesce(B.count,0) as count
  from (select distinct well_id from msm_energywellsphipsi where subspace_id = ss_id) A 
       left join 
       (select well_id1 as well_id, count(*) as count from MSM_Transitions group by well_id) B 
  on A.well_id = B.well_id;

  select well_id into min_well from TempWellStats order by count limit 1;
  select array_agg(centroid order by feature_id) into min_well_features
  from MSM_EnergyWellsPhiPsi
  where subspace_id = ss_id and well_id = min_well and ref_id = 1;

  -- for each expjob find the corresponding trjectory that just ran
  -- for each trajectory check whether the trajectory has reached one of the
  -- destination wells yet. If so, schedule a trajectory to start from an energy well.
  -- Otherwise, continue the current trajectory.

  for e_id in (select expjob_id from ExpJobs E where E.subspace_id = ss_id) loop
    select T.trj_id, M.well_id1, T.tran_id into trj_id, starting_well, tran_id
    from Trajectories T, MSM_Transitions M
    where     T.expjob_id = e_id 
          and T.tran_id   = M.tran_id
    order by T.trj_id desc limit 1;


    if trj_id is not Null then
      if exists (select 1 from MSM_Transitions where trj_id2 = trj_id limit 1) then
        insert into MSM_Transitions (subspace_id, expjob_id, well_id1, trj_id1) values (ss_id, e_id, min_well, trj_id);
        start_from_trj = Null;
        features = min_well_features;
      else
        start_from_trj = trj_id;
        features = Null;
      end if;
    else
      insert into MSM_Transitions (subspace_id, expjob_id, well_id1) values (ss_id, e_id, min_well);
      start_from_trj = Null;
      features = min_well_features;
    end if;

    select M.tran_id into tran_id 
    from MSM_Transitions M 
    where subspace_id = ss_id and expjob_id = e_id 
    order by M.tran_id desc limit 1;




    perform queue_gearman_job(e_id, features, 'straight_md', null, iter, start_from_trj, tran_id);
    update ExpJobs
    set blocked = false
    where expjob_id = e_id;
    num_jobs = num_jobs + 1;
  end loop;

  return num_jobs;
end
$$ language plpgsql;

--
--drop function if exists schedule_straight_md_jobs(ss_id int);
--create or replace function schedule_straight_md_jobs(ss_id int)
--returns int as $$
--declare
--  iter  int;
--  a_id  int;
--  st_id int;
--  c     int;
--  greediness float;
--begin
--
--  select get_state(ss_id) into st_id;
--  --a_id = 18;
--
--  select acc_score into greediness
--  from RL_History
--  where subspace_id = ss_id and
--        acc_score is not null
--  order by iter_id desc
--  limit 1;
--
--
--  select get_action(ss_id, st_id, greediness) into a_id;
--
--  insert into RL_History (subspace_id, state_id, action_id) values (ss_id, st_id, a_id);
--
--  select iter_id into iter
--  from RL_History
--  where subspace_id = ss_id and
--        acc_score is null
--  order by iter_id desc
--  limit 1;
--
--  select take_action(ss_id, a_id, iter) into c;
--
--  return c;
--end
--$$ language plpgsql;
--


drop function if exists take_action(ss_id bigint, a_id bigint, iter int);
create or replace function take_action(ss_id bigint, a_id bigint, iter int)
returns int as $$
declare
  r_id int[];
  rpd int[];
  nd int;
  c int = 0;
  e_id bigint;
  peak_ratio float = 0;
  nr int;
  o int;
  f float[];
  prev_iter int;
  start_from_trj int;
begin
  perform generate_regions_to_sample(iter, ss_id);

  select A.peak_ratio into peak_ratio from rl_actions A where A.action_id = a_id;

  for e_id in (select expjob_id from ExpJobs E where E.subspace_id = ss_id) loop

    if random() < peak_ratio then
      select count(*) into nr from Peaks_to_sample;
      select get_random_number(1,nr/4) into o;
    
      select region_id into r_id
      from Peaks_to_sample
      order by score limit 1 offset o;
    else
      select count(*) into nr from Wells_to_sample;
      select get_random_number(1,nr/4) into o;
 
      select region_id into r_id
      from Wells_to_sample
      order by score limit 1 offset o;
    end if;
 
    select count(*) into nd from unnest(r_id);
  
    select regions_per_dim into rpd
    from HistogramSpec
    where num_dims = nd;
 
    select array_agg(regionid_to_rad(r_id[A.a], rpd[A.a]) order by A.a)
    into f
    from generate_series(1,nd) A;
  
    start_from_trj = null;

    select T.trj_id into start_from_trj from Trajectories T
    where expjob_id = e_id
    order by T.trj_id desc limit 1;

    perform queue_gearman_job(e_id, f, 'straight_md', null, iter, start_from_trj);
    update ExpJobs
    set blocked = false
    where expjob_id = e_id;
    c = c+1;
  end loop;

  return c;
end
$$ language plpgsql;

drop function if exists assess_action(bigint);
create or replace function assess_action(ss_id bigint)
returns void  as $$
declare
  st_id int;
  a_id int;
  prev_score float;
  prev_st_id bigint;
  prev_action_id bigint;
  reward     float;
  delta_q    float;
  alpha      float = 1;
  gamma      float = 1;
  prev_q     float;
  curr_q     float;
  iter      int;
  b  float;
  num_states int = 0;
  cos_sim    float;
begin

  select count(*) into num_states
  from RL_History
  where subspace_id = ss_id;

  if num_states > 0 then
    select state_id,   action_id,      coalesce(acc_score,0)
    into   prev_st_id, prev_action_id, prev_score
    from RL_History
    where subspace_id = ss_id and
          acc_score is not null
    order by iter_id desc
    limit 1;

    select iter_id,state_id,action_id
    into   iter,st_id, a_id
    from RL_History
    where subspace_id = ss_id and
          acc_score is null
    order by iter_id desc
    limit 1;

    perform update_histograms(iter);

    select compute_accuracy(iter) into cos_sim;
    select compute_bumpiness(iter) into b;

    update RL_History 
    set acc_score = cos_sim,
        bumpiness = b
    where iter_id = iter;

    if prev_score is not null then
    
      select coalesce(q_value,0) into prev_q
      from RL_QTable
      where subspace_id = ss_id and 
            action_id = prev_action_id and 
            state_id = prev_st_id;
    
      select coalesce(q_value,0) into curr_q
      from RL_QTable
      where subspace_id = ss_id and 
            action_id   = a_id and 
            state_id    = st_id;
    
   
      reward = cos_sim - prev_score;

      delta_q = alpha * (reward + gamma * prev_q - curr_q);
    
      if delta_q is null then
        delta_q = 0;
      end if;

      update RL_QTable
      set q_value = q_value + delta_q
      where subspace_id = ss_id and 
            action_id   = a_id and 
            state_id    = st_id;
    end if;
  end if;
end
$$ language plpgsql;



drop function if exists regionid_to_rad(r_id int, num_regions int);
create or replace function regionid_to_rad(r_id int, num_regions int)
returns float as $$
begin
  return (r_id + random()) * (2*pi()/num_regions) - pi();
end
$$ language plpgsql;

drop function if exists resample_blocking(gmd_host text);
create or replace function resample_blocking(gmd_host text)
returns int as $$
declare
  num_trjs_consumed int;
  bulk_load_size    int = 4;
  e_id              bigint;
  ss_id             int;
  t_id              bigint;
  blocking          boolean;
  all_blocked       boolean;
  state_id          int;
  action_id         int;
  biased            boolean = false;
  done_waiting      boolean;
  num_jobs          int;
  inactive          boolean;
  num_trjs          int = 0;
  curr_trj          int;
begin

  select trj_id into t_id from FinishedTrajectories order by trj_id limit 1;
  if t_id is not null
  then

    select T.expjob_id,T.biased, E.subspace_id into e_id, biased, ss_id
    from Trajectories T,
         ExpJobs E
    where T.trj_id = t_id and
          T.expjob_id = E.expjob_id;

    if e_id is null then
      delete from FinishedTrajectories
      where trj_id = t_id;
      return 0;
    end if;



    update Subspaces
    set num_waits = num_waits - 1
    where subspace_id = ss_id;

    if not biased then
      update ExpJobs
      set num_iterations = num_iterations - 1,
          active   = (num_iterations - 1) > 0,
          blocked  = true
      where expjob_id = e_id;
    end if;
 
    select bool_and(blocked) into done_waiting
    from ExpJobs
    where subspace_id = ss_id;

    if done_waiting then
      select bool_and(not active) into inactive from ExpJobs where subspace_id = ss_id;

      --perform assess_action(ss_id);

      perform detect_transitions(ss_id);
      if inactive then
        delete from FinishedTrajectories
        where trj_id = t_id;

        return 0;
      end if;
      select schedule_straight_md_jobs(ss_id) into num_jobs;
      update Subspaces
      set num_waits = num_jobs,
          status    = 'straight_md'
      where subspace_id = ss_id;
    end if;

    delete from FinishedTrajectories
    where trj_id = t_id;
    select count(*) into num_trjs from  FinishedTrajectories;
  end if;



  return num_trjs;
end
$$ language plpgsql;



drop function if exists load_trajectories();
create or replace function load_trajectories()
returns int as $$
declare
  num_trjs_consumed int;
  bulk_load_size int = 1;
  num_new_jobs int = 0;
  rec record;
  e_id bigint;
begin

  if exists (select * from FinishedTrajectories limit 1) 
  then
    select expjob_id into e_id
    from FinishedTrajectories order by trj_id limit 1;

    update ExpJobs 
    set num_iterations = num_iterations - 1,
        active   = (num_iterations - 1) >= 0
    where expjob_id = e_id;

    create temp table tables_to_load as
    select *
    from FinishedTrajectories 
    where expjob_id = e_id 
    order by trj_id limit bulk_load_size;

    select count(*) into num_trjs_consumed from tables_to_load;

    perform process_new_data(e_id);

    delete from FinishedTrajectories
    where hs_t_name in (select hs_t_name from tables_to_load);
  end if;
  return num_trjs_consumed;
end 
$$ language plpgsql;


drop function if exists get_state(ss_id int);
create or replace function get_state(ss_id int)
returns int as $$
declare
  st_id int;
  b float;
begin
  select bumpiness into b from RL_History where subspace_id = ss_id and bumpiness is not null
  order by iter_id desc limit 1;
  st_id = floor(b*10);
  if not exists (select 1 from RL_History where subspace_id = ss_id and 
                                                state_id    = st_id limit 1) and st_id is not null
  then
    insert into RL_QTable (subspace_id, state_id, action_id, q_value)
    select ss_id,st_id,A.action_id,0 from RL_Actions A;

  end if;

  return st_id;
end 
$$ language plpgsql;


--
--drop function if exists get_state(ss_id int);
--create or replace function get_state(ss_id int)
--returns int as $$
--declare
--  vec_size float;
--  next_state_id int;
--  nearest_state_id int;
--  nearest_state_dist float;
--  threshold float default 0.100;
--begin
--  create temp table HistogramTemp as
--  select region_id,sum(cast(n as float)) as n
--  from histograms H, trajectories T, expjobs E 
--  where H.trj_id = T.trj_id and T.expjob_id = E.expjob_id and E.subspace_id = ss_id
--  group by region_id;
--
--  select sqrt(sum(n*n)) into vec_size from HistogramTemp;
--
--  update HistogramTemp set n = n/vec_size;
--
--  select state_id, 1-sum(S.n*T.n) as dist into nearest_state_id, nearest_state_dist
--  from RL_States S, HistogramTemp T
--  where S.region_id = T.region_id
--  group by state_id
--  order by dist
--  limit 1;
--
--  if (nearest_state_dist is null) or (nearest_state_dist > threshold) then
--
--    select coalesce(max(state_id),0) + 1 into next_state_id from RL_States;
--
--    insert into RL_States
--    select next_state_id as state_id, ss_id, H.region_id, H.n from HistogramTemp H;
--
--    insert into RL_QTable
--    select ss_id,next_state_id,A.action_id,random() from RL_Actions A;
--    return next_state_id;
--  end if;
--
--  return nearest_state_id;
--end
--$$ language plpgsql;
--
