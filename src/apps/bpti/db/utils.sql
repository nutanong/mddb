drop function if exists jsonify(text);

drop function if exists generate_job_json(bigint, bigint, text, float[]);
drop function if exists tor_distl1(float[], float[], float) cascade;
drop function if exists vec_length(float8, float8, float8);
drop function if exists vec_dot(float8, float8, float8, float8, float8, float8);
drop function if exists vec_cross(in x1 float8, in y1 float8, in z1 float8,
                          in x2 float8, in y2 float8, in z2 float8,
                          out x float8, out y float8, out z float8);

drop function if exists vector_angle(in x1 float8, in y1 float8, in z1 float8,
                                     in x2 float8, in y2 float8, in z2 float8,
                                     out angle float8);

drop function if exists dihedral_angle(in x1 float8, in y1 float8, in z1 float8,
                               in x2 float8, in y2 float8, in z2 float8,
                               in x3 float8, in y3 float8, in z3 float8,
                               in x4 float8, in y4 float8, in z4 float8,
                               out angle float8);
drop function if exists get_random_number(integer, integer);

drop function if exists update_Trajectories(bigint, bigint);
drop function if exists queue_gearman_job(bigint, bigint, bigint, float[], text, text);


drop function if exists generate_state_crd(numstates int);
drop function if exists incr_crd_file(int, int, int, float, float, float);
drop type if exists crd_file_state cascade;

drop function if exists generate_state_pdb(numstates int);
drop function if exists incr_pdb_file(bigint, int, int, float, float, float);
drop type if exists pdb_file_state cascade;


drop view if exists LockView;
create view LockView as
   select 
     pg_stat_activity.datname,pg_class.relname,pg_locks.transactionid, pg_locks.mode, pg_locks.granted,
     pg_stat_activity.usename,substr(pg_stat_activity.current_query,1,30), pg_stat_activity.query_start, 
     age(now(),pg_stat_activity.query_start) as "age", pg_stat_activity.procpid 
   from pg_stat_activity,pg_locks left 
     outer join pg_class on (pg_locks.relation = pg_class.oid)  
   where pg_locks.pid=pg_stat_activity.procpid and
         not granted
   order by query_start;

create or replace function id2phipsi(f float[]) returns float[] as $$
declare
  res float[];
  nd int;
begin
  select count(*) into nd from unnest(f);

  select array_agg(B.bucket_start + random()*(B.bucket_end - B.bucket_start) order by feature_id) into res
  from Buckets B 
  where num_dims = nd and bucket_id = f[feature_id];

  return res;
end
$$ language plpgsql;


create or replace function jsonify(st_in text) returns text as $$
declare
  st_out text;
begin
  st_out := array_to_string(string_to_array(st_in, E'\n'), E'","');
  st_out := E'["' || st_out || E'"]';
  return st_out;
end
$$ language plpgsql;

create or replace function generate_job_json(
  t_id           bigint,
  e_id           bigint, 
  p_seq          text,
  feature        float[],
  start_from_trj int
)
returns text
as $$
declare
  line      text;
  lines     text[];
  out_text  text;
  rec       record;
  ns_sim    int;
  ts_pico   float;
  tsav_freq int;
  simu      text;
  sim_ff text;
  p_id   int;
begin
  lines = array_append(lines, E'"e_id": '         || e_id);
  lines = array_append(lines, E'"trj_id": '       || t_id);
  lines = array_append(lines, E'"protein_seq": "' || p_seq   || '"');

  select nsteps_sim, timestep_pico, trjsave_freq, simulator 
  into ns_sim, ts_pico, tsav_freq, simu
  from ExpJobs where expjob_id = e_id;

  lines = array_append(lines, E'"nstep_simulation": '  || ns_sim);
  lines = array_append(lines, E'"timestep_pico": '     || ts_pico);
  lines = array_append(lines, E'"trj_save_freq": '     || tsav_freq);
  lines = array_append(lines, E'"simulator": "'        || simu || '"');

  select ff, S.protein_id into sim_ff, p_id
  from SimConfigs C, Subspaces S, ExpJobs E
  where S.config_id = C.config_id and 
        S.subspace_id = E.subspace_id and 
        E.expjob_id = e_id;

  lines = array_append(lines, E'"sim_ff": "'           || sim_ff || '"');
  lines = array_append(lines, E'"protein_id": '        || p_id);
  lines = array_append(lines, E'"start_from_trj": '    || start_from_trj);

  if not feature is Null then
    line = E'"phi_psi_array": [' || array_to_string(feature, ', ') || ']';
  else
    line = E'"phi_psi_array": []';
  end if;
  lines = array_append(lines, line);

  --for rec in (select param_name, param_value
  --            from MDJobParams where trj_id = t_id)
  --loop
  --  line = E'"' || rec.param_name || E'": "' || rec.param_value || '"';
  --  lines = array_append(lines,line);
  --end loop;

  out_text = E'{\n' || array_to_string(lines, E',\n') || E'\n}';

  return out_text;
end
$$ language plpgsql;



create or replace function generate_job_json(
  t_id           bigint,
  ss_id          bigint,
  e_id           bigint,
  feature_id     int,
  p_seq          text,
  phi_psi        float
)
returns text
as $$
declare
  line      text;
  lines     text[];
  out_text  text;
  rec       record;
  ns_sim    int;
  ts_pico   float;
  tsav_freq int;
  simu      text;
  sim_ff    text;
  p_id      int;
begin


  lines = array_append(lines, E'"trj_id": '        || t_id);
  lines = array_append(lines, E'"ss_id": '         || ss_id);
  lines = array_append(lines, E'"e_id": '          || e_id);
  lines = array_append(lines, E'"feature_id": '    || feature_id);
  lines = array_append(lines, E'"protein_seq": "'  || p_seq || '"');
  lines = array_append(lines, E'"phi_psi_array": ' || phi_psi);
  

  select nsteps_sim, timestep_pico, trjsave_freq, simulator
  into ns_sim, ts_pico, tsav_freq, simu
  from ExpJobs where expjob_id = e_id;

  lines = array_append(lines, E'"nstep_simulation": '  || ns_sim);
  lines = array_append(lines, E'"timestep_pico": '     || ts_pico);
  lines = array_append(lines, E'"trj_save_freq": '     || tsav_freq);
  lines = array_append(lines, E'"simulator": "'        || simu || '"');

  select ff, S.protein_id into sim_ff, p_id
  from SimConfigs C, Subspaces S, ExpJobs E
  where S.config_id = C.config_id and 
        S.subspace_id = E.subspace_id and 
        E.expjob_id = e_id;

  lines = array_append(lines, E'"sim_ff": "'           || sim_ff || '"');

  lines = array_append(lines, E'"protein_id": '        || p_id);


  out_text = E'{\n' || array_to_string(lines, E',\n') || E'\n}';

  return out_text;
end
$$ language plpgsql;


create or replace function tor_distl1(v1 float[], v2 float[], diff_width float)
returns float as $$
declare
  l int;
  i int;
  d1 float;
  d2 float;
  dist float = 0;
begin
  l = 2;
  for i in (select * from generate_series(1,l)) loop
    d1 = abs(v1[i] - v2[i]);
    d2 = diff_width - d1;
    if d1 < d2 then
      dist = dist + d1;
    else
      dist = dist + d2;
    end if;
  end loop;

  return dist;
end;
$$ language plpgsql;








create function vec_length(x float8, y float8, z float8) returns float8
as $$
  begin
    return sqrt(x*x+y*y+z*z);
  end;
$$ language plpgsql;

create function vec_dot(x1 float8, y1 float8, z1 float8,
                        x2 float8, y2 float8, z2 float8) returns float8
as $$
  begin
    return x1*x2+y1*y2+z1*z2;
  end;
$$ language plpgsql;


create function vec_cross(in x1 float8, in y1 float8, in z1 float8,
                          in x2 float8, in y2 float8, in z2 float8,
                          out x float8, out y float8, out z float8)
as $$
  begin
    x := (y1*z2-z1*y2);
    y := (z1*x2-x1*z2);
    z := (x1*y2-y1*x2);
  end;
$$ language plpgsql;


create function vector_angle(in x1 float8, in y1 float8, in z1 float8,
                             in x2 float8, in y2 float8, in z2 float8,
                             out angle float8)
as $$
  begin
  angle := acos(vec_dot(x1,y1,z1,x2,y2,z2) /
               (vec_length(x1,y1,z1)*vec_length(x2,y2,z2)));
  end
$$ language plpgsql;



create function dihedral_angle(in x1 float8, in y1 float8, in z1 float8,
                               in x2 float8, in y2 float8, in z2 float8,
                               in x3 float8, in y3 float8, in z3 float8,
                               in x4 float8, in y4 float8, in z4 float8,
                               out angle float8)
as $$
  declare
    v1_x float8; v1_y float8; v1_z float8;
    v2_x float8; v2_y float8; v2_z float8;
    v3_x float8; v3_y float8; v3_z float8;
    n1_x float8; n1_y float8; n1_z float8;
    n2_x float8; n2_y float8; n2_z float8;
  begin

  v1_x = x2-x1;
  v1_y = y2-y1;
  v1_z = z2-z1;

  v2_x = x3-x2;
  v2_y = y3-y2;
  v2_z = z3-z2;

  v3_x = x4-x3;
  v3_y = y4-y3;
  v3_z = z4-z3;

  select (N1.v).x, (N1.v).y, (N1.v).z,
         (N2.v).x, (N2.v).y, (N2.v).z
    into n1_x, n1_y, n1_z, n2_x, n2_y, n2_z
  from (select vec_cross(v1_x, v1_y, v1_z,
                         v2_x, v2_y, v2_z) as v) as N1,
       (select vec_cross(v2_x, v2_y, v2_z,
                         v3_x, v3_y, v3_z) as v) as N2;

  angle = atan2(vec_length(v2_x, v2_y, v2_z)*vec_dot(v1_x,v1_y,v1_z,n2_x,n2_y,n2_z),
                vec_dot(n1_x, n1_y, n1_z, n2_x, n2_y, n2_z));

  end
$$ language plpgsql;


create or replace function get_random_number(start_int integer, end_int integer) 
returns integer as $$
begin
    return trunc(random() * (end_int-start_int) + start_int);
end;
$$ LANGUAGE plpgsql;


create or replace function queue_gearman_job(
  e_id              bigint,
  feature_array     float[],
  task_id           text,
  biased_feature    int,
  iter              int,
  start_from_trj    int,
  tran_id           int
)
returns void
as $$
declare
  p_seq          text;
  trj_id         bigint;
  job_data       text;
  job_results    text;
  ss_id          bigint;
  gh             text;
  gp             int;
begin

  select trajectories_insert(e_id, biased_feature, iter, start_from_trj, tran_id) into trj_id;


  select subspace_id into ss_id
  from ExpJobs E
  where E.expjob_id = e_id;

  select protein_cname into p_seq 
  from Subspaces S, Proteins P 
  where S.subspace_id = ss_id and 
        P.protein_id = S.protein_id;

  select gmdhost, gmdport into gh, gp
  from Subspaces S
  where S.subspace_id = ss_id;

  job_data = generate_job_json(trj_id, e_id, p_seq, feature_array, start_from_trj);

  perform gman_servers_set(gh || ':' || gp);
  job_results = gman_do_background(task_id, job_data);
end
$$ language plpgsql;



create or replace function queue_gearman_job(
  task_id           text,
  job_data          text
)
returns void
as $$
declare
  p_seq          text;
  trj_id         bigint;
  job_results    text;
  ss_id          bigint;
  gh             text;
  gp             int;
begin

  select param_value from controlparams into gh where param_name = 'gmdhost';
  select param_value from controlparams into gp where param_name = 'gmdport';

  perform gman_servers_set(gh || ':' || gp);
  job_results = gman_do_background(task_id, job_data);
end
$$ language plpgsql;





create or replace function incr_crd_file(r_trj_id int, r_t int, r_atom_id int,
                                         x float, y float, z float)
returns text
as $$
declare
  tmp text;
begin
  select lpad(cast(atom_id as text), 5) ||
         lpad(cast(A.residue_id as text), 5) || E' ' ||
         rpad(A.residue_name,4) || E' '||
         rpad(A.atom_name, 4) ||
         lpad(trim(to_char(x, '9990.99999')), 9) ||
         lpad(trim(to_char(y, '9990.99999')), 9) ||
         lpad(trim(to_char(z, '9990.99999')), 9) || E' ' ||
         rpad(A.segment_name,4) || E' ' ||
         rpad(trim(to_char(A.residue_id,'9999')),4) ||
         lpad(trim(to_char(0.0, '99990.99999')), 10)
         into tmp
  from  AtomMetaTemp A
  where A.atom_id = r_atom_id;
  return tmp;
end
$$ language plpgsql
immutable;

create type crd_file_state as (
  atom_count  int,
  atom_coords text
);

create or replace function compose_crd_file(s crd_file_state) returns text
as $$
declare
  crd_template text;
begin
  crd_template :=
    E'* CRD FILE FOR ALANINE DIPEPTIDE\n' ||
    E'*  DATE:   ' || now() ||
    E' CREATED BY USER: ' || current_user ||E'\n' ||
    E'*\n' ||
    lpad(trim(to_char(s.atom_count, '99999')), 5) || E'\n' ||
    s.atom_coords || E'\n';

  return crd_template;
end
$$ language plpgsql;


create or replace function generate_state_crd(numstates int)
  returns text[] as $$
declare
  rowdata record;
  point record;
  temp text;
  temp1 crd_file_state;
  results text[];
  num_confs_avaliable int;
begin
  for point in (select * from Conf_to_sample order by score limit numstates)
  loop

    temp1.atom_count = 0;
    temp1.atom_coords = '';
    for rowdata in (select P.trj_id, P.t, P.atom_id, P.x, P.y, P.z from
                      AtomPositions P
                      where point.trj_id = P.trj_id and point.t = P.t
                      order by trj_id, t, atom_id)
    loop
      select incr_crd_file(rowdata.trj_id, rowdata.t,
                           rowdata.atom_id, rowdata.x, rowdata.y, rowdata.z) into temp;

      temp1 := row(temp1.atom_count+1, temp1.atom_coords || temp || E'\n');
    end loop;

    results := array_append(results, compose_crd_file(temp1));
  end loop;

  return results;
end
$$ language plpgsql;

-- Generate PDB

create type pdb_file_state as (
  atom_count      int,
  st              text,
  last_residue    text,
  last_residue_id int,
  trj_id          int,
  t               int
);

create or replace function compose_pdb_entry(atom_id int, atom_name text, residue_name text, 
                                             residue_id int, x float, y float, z float,
                                             segment_name text)
returns text
as $$
declare
  res text;
begin
  res = 
         'HETATM  ' ||
         lpad(cast(atom_id as text), 3) || E' '||
         rpad(atom_name, 4) || E' '||
         rpad(residue_name, 4) || E' '||
         lpad(cast(residue_id as text), 4) || E' ' ||
         lpad(trim(to_char(x, '990.999')), 11) || E' ' ||
         lpad(trim(to_char(y, '990.999')), 7) || E' ' ||
         lpad(trim(to_char(z, '990.999')), 7) || E' ' ||
         '   1.00  0.00 ' ||
         lpad(segment_name, 5) || E'  ';
  return res;
end;
$$ language plpgsql;

create or replace function incr_pdb_file(r_trj_id bigint, r_t int, r_atom_id int,
                                         x float, y float, z float)
returns pdb_file_state
as $$
declare
  ret pdb_file_state;
  rec record;
begin
  select
         'HETATM  ' ||
         lpad(cast(atom_id as text), 3) || E' '||
         rpad(A.atom_name, 4) || E' '||
         rpad(A.residue_name, 4) || E' '||
         lpad(cast(A.residue_id as text), 4) || E' ' ||
         lpad(trim(to_char(x, '990.999')), 11) || E' ' ||
         lpad(trim(to_char(y, '990.999')), 7) || E' ' ||
         lpad(trim(to_char(z, '990.999')), 7) || E' ' ||
         '   1.00  0.00 ' ||
         lpad(A.segment_name, 5) || E'  ' as st,
         A.residue_name as last_residue,
         A.residue_id as last_residue_id
         into rec
  from  AtomMetaTemp A
  where   A.atom_id = r_atom_id;

  ret.atom_count = 0;
  ret.st = rec.st;
  ret.last_residue_id = rec.last_residue_id;
  ret.last_residue = rec.last_residue;
  
  return ret;
end
$$ language plpgsql
immutable;

create or replace function compose_pdb_file(trj_id int,t int,st text) returns text
as $$
declare
  crd_template text;
begin
  crd_template :=
    -- Header
    E'REMARK trj_id = ' || trj_id || E'; t = ' || t || E'\n' ||
    E'REMARK  DATE:   ' || now() || '       ' ||
    E'CREATED BY USER: ' || current_user ||E'\n' ||
    -- Actual Content Coordinates
    st ||
    -- Terminator
    --E'TER        ' || (s.atom_count+1) || 
    --E'       ' || s.last_residue || 
    --E'       ' || s.last_residue_id || E'\n' ||
    E'END\n';

  return crd_template;
end
$$ language plpgsql;



create or replace function generate_pdb(trj_id int, ts int)
  returns text as $$
declare
  st      text;
begin
  select 
         E'REMARK trj_id = ' || trj_id || E'; t = ' || ts || E'\n' ||
         E'REMARK  DATE:   ' || now() || '       ' ||
         E'CREATED BY USER: ' || current_user ||E'\n' ||
         string_agg(compose_pdb_entry(P.atom_id, M.atom_name, M.residue_name, M.residue_id, 
                                      P.x, P.y, P.z, M.segment_name), E'\n' order by P.atom_id) ||
         E'END\n'
  into st
  from AtomPositions P, Proteinatoms M, Subspaces S, ExpJobs E, Trajectories T
  where S.subspace_id = E.subspace_id and
        E.expjob_id = T.expjob_id and
        T.trj_id = trj_id and
        P.trj_id = trj_id and 
        P.t = ts and
        M.atom_id = P.atom_id and
        M.config_id  = S.config_id and
        M.protein_id = S.protein_id;

  return st;
end
$$ language plpgsql;

create or replace function generate_pdb(conf_id int)
  returns text as $$
declare
  t int;
  trj_id int;
  
begin
  select A.t, A.trj_id into t, trj_id from HistoBinAssignments A where A.conf_id = conf_id;
  return generate_pdb(trj_id,t);
end
$$ language plpgsql;




create or replace function find_gpu(gwn text) returns text
as $$
declare
  id bigint;
begin
  lock table GPUs in exclusive mode;

  select device_id into id from GPUs where occupied = false order by device_id limit 1;

  if id is not null then
    update GPUs set (gw_name, occupied) = (gwn,true) where device_id = id;
    return id;
  else
    return -1;
  end if;
end
$$ language plpgsql;



drop function if exists generate_histogram(s_id int, duration_hour float, x int, y int);
create or replace function generate_histogram(s_id int, duration_hour float, x int, y int)
returns void as $$
begin
  perform generate_histogram(s_id, duration_hour, x, y, '/tmp/');
end
$$ language plpgsql;

drop function if exists generate_histogram(s_id int, duration_hour float, x int, y int, z int);
create or replace function generate_histogram(s_id int, duration_hour float, x int, y int, z int)
returns void as $$
begin
  perform generate_histogram(s_id, duration_hour, x, y, z, '/tmp/');
end
$$ language plpgsql;


drop function if exists combine_histograms(s_id int, duration_hour float);
create or replace function combine_histograms(s_id int, duration_hour float)
returns void as $$
declare
  t_start timestamp with time zone;
begin
  select min(t1) into t_start 
  from Benchmark B, Trajectories T, ExpJobs E
  where B.trj_id = T.trj_id and
        T.expjob_id = E.expjob_id and
        E.subspace_id = s_id;

  drop table if exists CombinedHistogram;
  create temp table CombinedHistogram as
  select id_array, region_id,
         sum(n) as n,
         sum(s) as s
  from Histograms H, Trajectories T, Benchmark B, ExpJobs E
  where T.expjob_id = E.expjob_id and
        E.subspace_id = s_id and
        H.trj_id = T.trj_id and 
        B.trj_id = T.trj_id and
        extract(epoch from t3-t_start) < duration_hour*3600.0
  group by id_array, region_id;
end
$$ language plpgsql;



drop function if exists write_histogram2d(x int, y int, outprefix text);
create or replace function write_histogram2d(x int, y int, outprefix text)
returns void as $$
begin
  drop table if exists ProjectedHistogram;
  create temp table ProjectedHistogram as
  select id_array[x] as a, id_array[y] as b, sum(n) as n, sum(s) as s
  from CombinedHistogram
  group by id_array[x], id_array[y];

  execute 'copy ProjectedHistogram to ' || 
          quote_literal(outprefix || 'histogram.csv') || 
          ' with delimiter ' || quote_literal(',');

  drop table if exists ProjectedHistogramRegion;
  create temp table ProjectedHistogramRegion as
  select region_id[x] as a, region_id[y] as b, sum(n) as n, sum(s) as s, sum(n)*sum(s) as score
  from CombinedHistogram
  group by region_id[x], region_id[y];

  execute 'copy ProjectedHistogramRegion to ' ||
          quote_literal(outprefix || 'histo_scores.csv') ||
          ' with delimiter  ' || quote_literal(',');
end
$$ language plpgsql;

drop function if exists write_histogram3d(x int, y int, z int, outprefix text);
create or replace function write_histogram3d(x int, y int, z int, outprefix text)
returns void as $$
begin
  drop table if exists ProjectedHistogram;
  create temp table ProjectedHistogram as
  select id_array[x] as a, id_array[y] as b, id_array[z] as c, sum(n) as n, sum(s) as s
  from CombinedHistogram
  group by id_array[x], id_array[y], id_array[z];

  execute 'copy ProjectedHistogram to ' ||
          quote_literal(outprefix || 'histogram.csv') ||
          ' with delimiter ' || quote_literal(',');

  drop table if exists ProjectedHistogramRegion;
  create temp table ProjectedHistogramRegion as
  select region_id[x] as a, region_id[y] as b, region_id[z] as c, sum(n) as n, sum(s) as s, sum(n)*sum(s) as score
  from CombinedHistogram
  group by region_id[x], region_id[y], region_id[z];

  execute 'copy ProjectedHistogramRegion to ' ||
          quote_literal(outprefix || 'histo_scores.csv') ||
          ' with delimiter  ' || quote_literal(',');
end
$$ language plpgsql;



drop function if exists generate_histogram(s_id int, duration_hour float, x int, y int, outprefix text);
create or replace function generate_histogram(s_id int, duration_hour float, x int, y int, outprefix text)
returns void as $$
begin
  perform combine_histograms(s_id, duration_hour);
  perform write_histogram2d(x,y,outprefix);
end
$$ language plpgsql;



drop function if exists generate_histogram(s_id int, duration_hour float, x int, y int, z int, outprefix text);
create or replace function generate_histogram(s_id int, duration_hour float, x int, y int, z int, outprefix text)
returns void as $$
begin
  perform combine_histograms(s_id, duration_hour);
  perform write_histogram3d(x,y,z,outprefix);
end
$$ language plpgsql;



drop function if exists generate_1d_histograms(ss_id int, resolution int);
create or replace function generate_1d_histograms(ss_id int, resolution int)
returns void as $$
declare
  iter int;
  res int;
  bins_per_region   int;

begin
  select cast(param_value as int) into res
  from ControlParams where param_name = 'histo_res_lo';

  bins_per_region = res/resolution;


  select iter_id into iter
  from RL_History
  where subspace_id = ss_id and
        acc_score is not null
  order by iter_id desc
  limit 1;

  drop table if exists TempHistograms;
  create table TempHistograms as
  select F.feature_id,bucket_id[F.feature_id]/bins_per_region as key,sum(n) as n,sum(s) as s
  from IterHistograms H, Subspaces S, ProteinFeatures F
  where H.iter_id = iter and
        S.subspace_id = ss_id and
        F.protein_id  = S.protein_id
  group by F.feature_id, key;
  
end
$$ language plpgsql;

create or replace function generate_2d_histogram(ss_id int, x int, y int, resolution int)
returns void as $$
declare
  iter int;
  res int;
  bins_per_region   int;

begin
  select cast(param_value as int) into res
  from ControlParams where param_name = 'histo_res_hi';

  bins_per_region = res/resolution;

  select iter_id into iter
  from RL_History
  where subspace_id = ss_id and
        acc_score is not null
  order by iter_id desc
  limit 1;

  drop table if exists TempHistogram;
  create table TempHistogram as
  select bucket_id[x]/bins_per_region as x_key,
         bucket_id[y]/bins_per_region as y_key,
         sum(n) as n,sum(s) as s
  from IterHistogramsHI H
  where H.iter_id = iter 
  group by x_key, y_key;
  
end
$$ language plpgsql;

create or replace function generate_2d_histogram(ss_id int, st_id int, x int, y int, resolution int)
returns void as $$
declare
  iter int;
  histogram_max_res int;
  bins_per_region   int;

begin
  select cast(param_value as int) into histogram_max_res
  from ControlParams where param_name = 'histogram_max_res';

  bins_per_region = histogram_max_res/resolution;

  select iter_id into iter
  from RL_History
  where subspace_id = ss_id and
        acc_score is not null and
        state_id = st_id
  order by iter_id desc
  limit 1;

  drop table if exists TempHistogram;
  create table TempHistogram as
  select bucket_id[x]/bins_per_region as x_key,
         bucket_id[y]/bins_per_region as y_key,
         sum(n) as n,sum(s) as s
  from IterHistogramsHI H
  where H.iter_id = iter 
  group by x_key, y_key;
  
end
$$ language plpgsql;


create or replace function cache_load(iter bigint, acc_type text, acc_note text)
returns void as $$
declare
  evt_iter int;
  n int;
  ss_id int;
begin

  if iter is null then
    return;
  end if;

  if exists (select 1 from CacheTimestamps where iter_id = iter) then
    update CacheTimestamps 
    set ts=statement_timestamp(), 
        access_type = acc_type,
        note = acc_note
    where iter_id = iter;
    insert into CacheLogs (iter_id, ts, access_type, note,hit) 
    values (iter,statement_timestamp(), acc_type, acc_note, true);

    return;
  end if;

  select subspace_id into ss_id from RL_History R where R.iter_id = iter;

  perform cache_evict(ss_id);

  insert into IterHistogramsHICache
  select * from IterHistogramsHI where iter_id = iter;


  insert into IterHistogramsCache
  select * from IterHistograms where iter_id = iter;

  insert into PartitionedHistogramsCache
  select * from PartitionedHistograms where iter_id = iter;

  insert into CacheTimestamps (iter_id, ts, access_type, note) values (iter,statement_timestamp(), acc_type, acc_note);
  insert into CacheLogs (iter_id, ts, access_type, note,hit) 
  values (iter,statement_timestamp(), acc_type, acc_note,false);
end
$$ language plpgsql;

create or replace function cache_load_from_temps(iter bigint, acc_type text, acc_note text)
returns void as $$
declare
  evt_iter int;
  n int;
  ss_id int;
begin

  if iter is null then
    return;
  end if;

  if exists (select 1 from CacheTimestamps where iter_id = iter) then
    update CacheTimestamps
    set ts=statement_timestamp(),
        access_type = acc_type,
        note = acc_note
    where iter_id = iter;
    insert into CacheLogs (iter_id, ts, access_type, note,hit)
    values (iter,statement_timestamp(), acc_type, acc_note, true);

    return;
  end if;

  select subspace_id into ss_id from RL_History R where R.iter_id = iter;

  perform cache_evict(ss_id);

  insert into IterHistogramsHICache
  select * from TempIterHistogramsHI;


  insert into IterHistogramsCache
  select * from TempIterHistograms;

  insert into PartitionedHistogramsCache
  select * from TempPartitionedHistograms;


  insert into CacheTimestamps (iter_id, ts, access_type, note) values (iter,statement_timestamp(), acc_type, acc_note);
  insert into CacheLogs (iter_id, ts, access_type, note,hit)
  values (iter,statement_timestamp(), acc_type, acc_note,false);
end
$$ language plpgsql;



create or replace function cache_evict(ss_id int)
returns void as $$
declare
  cache_size int;
  e_iter int;
  n int;
begin
  select cast(param_value as int) into cache_size
  from ControlParams where param_name = 'cache_size';

  select count(*) into n 
  from CacheTimestamps C, rl_history R 
  where C.iter_id = R.iter_id;-- and R.subspace_id = ss_id;

  while n >= cache_size loop
    select C.iter_id into e_iter 
    from CacheTimestamps C, rl_history R 
    where C.iter_id = R.iter_id;-- and R.subspace_id = ss_id;

    delete from IterHistogramsCache where iter_id = e_iter;
    delete from IterHistogramsHICache where iter_id = e_iter;
    delete from PartitionedHistogramsCache where iter_id = e_iter;
    delete from CacheTimestamps where iter_id = e_iter;

    select count(*) into n 
    from CacheTimestamps C, rl_history R 
    where C.iter_id = R.iter_id;-- and R.subspace_id = ss_id;

  end loop;
end
$$ language plpgsql;



