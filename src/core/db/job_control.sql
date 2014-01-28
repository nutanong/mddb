drop table if exists  Workers;
drop table if exists  Benchmark cascade;
drop table if exists  ClusterConfig;
drop table if exists  MachineInfo;
drop table if exists  JobQueue;
drop table if exists  LoadingQueue;
drop view if  exists  JobRunningRecords;
drop table if exists  JobExecutionRecords;
drop table if exists  JobCompletionRecords;
drop table if exists  JobEntries;

create table Workers (
  worker_id     serial,
  worker_host   text,
  worker_name   text,
  process_id    int,
  resource_name text,
  primary key   (worker_id)
);


create table MachineInfo (
  host_name   text,
  key         text,
  value       text
);


drop function if exists workers_insert(text);
create or replace function workers_insert(h text)
returns int as $$
declare
  res int;
begin
  lock table Workers in exclusive mode;
  insert into Workers (worker_host) values (h);
  select max(worker_id) into res from Workers;
  return res;

end;
$$ language plpgsql;


create table JobEntries (
  jq_entry_id  bigint,
  job_id       int,
  primary key  (jq_entry_id)
);

create table JobQueue (
  jq_entry_id  serial,
  data         text,
  t_queue      timestamp with time zone,
  primary key  (jq_entry_id)
);

create table LoadingQueue (
  jq_entry_id   bigint,
  data          text,
  primary key   (jq_entry_id)
);

create table JobExecutionRecords (
  jq_entry_id    bigint,
  worker_id      int,
  session_dir    text,
  deployment_id  int,
  run_id         int,
  t_execution    timestamp with time zone,
  primary key    (jq_entry_id)
);

create table JobCompletionRecords (
  jq_entry_id    bigint,
  alright        bool,
  path           text,
  t_completion   timestamp with time zone,
  primary key    (jq_entry_id)
);

create view JobRunningRecords as
select E.jq_entry_id, E.worker_id, E.session_dir, E.deployment_id, E.run_id
from JobExecutionRecords E left join JobCompletionRecords C
on E.jq_entry_id = C.jq_entry_id
where C.jq_entry_id is null;


drop function if exists jobqueue_insert(text);
create or replace function jobqueue_insert(json_st text)
returns int as $$
declare
  max_id    int;
begin
  select jobqueue_insert(null, json_st) into max_id;

  return max_id;
end;
$$ language plpgsql;

drop function if exists jobqueue_insert(bigint, text);
create or replace function jobqueue_insert(j_id bigint, json_st text)
returns int as $$
declare
  max_id    int;
begin
  lock table JobQueue in exclusive mode;
  insert into JobQueue (data, t_queue)
  values (json_st, statement_timestamp());

  select max(jq_entry_id) into max_id from JobQueue;

  insert into JobEntries (jq_entry_id, job_id)
  values (max_id, j_id);
  return max_id;
end;
$$ language plpgsql;

drop function if exists jobqueue_dequeue(n int);
create or replace function jobqueue_dequeue(n int)
returns text[] as $$
declare
  ret text[];
begin
  lock table JobQueue in exclusive mode;

  create temp table SelectedJobs as
  select jq_entry_id, data
  from JobQueue
  order by jq_entry_id
  limit n;

  delete from JobQueue
  using SelectedJobs 
  where JobQueue.jq_entry_id = SelectedJobs.jq_entry_id;
  
  select array_agg(compose_jq_entry(jq_entry_id,data) order by jq_entry_id) into ret from SelectedJobs;

  return ret;
end;
$$ language plpgsql;

drop function if exists compose_jq_entry(jq_entry_id int, data text);
create or replace function compose_jq_entry(jq_entry_id int, data text)
returns text as $$
declare
  ret text;
begin
  ret = '{"jq_entry_id": ' || jq_entry_id || ', "data":' || data || '}';
  return ret;
end;
$$ language plpgsql;


drop function if exists loadingqueue_dequeue(n int);
create or replace function loadingqueue_dequeue(n int)
returns text[] as $$
declare
  ret text[];
begin
  lock table LoadingQueue in exclusive mode;

  create temp table SelectedJobs as
  select jq_entry_id, data
  from LoadingQueue
  order by jq_entry_id
  limit n;

  delete from LoadingQueue
  using SelectedJobs
  where LoadingQueue.jq_entry_id = SelectedJobs.jq_entry_id;

  select array_agg(data order by jq_entry_id) into ret from SelectedJobs;

  return ret;
end;
$$ language plpgsql;


