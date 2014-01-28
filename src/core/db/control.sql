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



drop function if exists mdq_pop();
create or replace function mdq_pop()
returns MDQEntry as $$
declare
  mdq_row        record;
  mdq_entry      MDQEntry;
  res            MDQEntry;
  term           text;
  s              text;
begin
  if not exists (select 1 from MDQueue limit 1) then
    return Null;
  end if;

  while True loop
    select Q.entry_id, Q.entry into mdq_row
    from MDQueue Q
    order by (Q.entry).num_terms_left, (Q.entry).expanded_terms limit 1;


    if found then
      delete from MDQueue where entry_id = mdq_row.entry_id;
      mdq_entry = mdq_row.entry;

      if mdq_entry.expanded_terms = Null then
        mdq_entry.expanded_terms = '{}'::text[];
      end if;
  
      if (mdq_row.entry).num_terms_left > 0 then
        mdq_entry.num_terms_left = (mdq_row.entry).num_terms_left - 1;
        term = (mdq_row.entry).rest_of_terms[1];
        mdq_entry.rest_of_terms = (mdq_row.entry).rest_of_terms[2:(mdq_row.entry).num_terms_left];
        if substring(term from 1 for 1) = '%' then
          for s in (select symbol from AminoAcids where property similar to term order by symbol)
          loop
            mdq_entry.expanded_terms = (mdq_row.entry).expanded_terms || s;
            insert into MDQueue (entry) values (mdq_entry);
          end loop;
        else
          mdq_entry.expanded_terms  = (mdq_row.entry).expanded_terms || term;
          insert into MDQueue (entry) values (mdq_entry);
        end if;
      end if;
    end if;
    exit when (mdq_row.entry).num_terms_left <= 0;
  end loop;
  res = mdq_row.entry;
  return res;
end;
$$ language plpgsql;


--drop function if exists mdq_pop();
--create or replace function mdq_pop()
--returns text as $$
--declare
--  heads text[];
--  tails text[];
--  term  text;
--  id    int;
--  l     int = 1;
--  ne    int;
--  ni    int;
--  ns    int;
--  ui    int;
--  p_id  int;
--  po    text;
--  sim   text;
--  c_id  int;
--  tl    text[];
--begin
--  while l > 0
--  loop
--
--    -- Dequeue: We order by num_terms_left so that we do depth first expansion
--    -- to save space.
--    select entry_id, protein_id, expanded_terms, rest_of_terms, num_terms_left, 
--           num_expjobs, num_iterations, num_steps_sim, user_id, policy,
--           simulator, config_id, trj_list
--    into id, p_id, heads, tails, l, ne, ni, ns, ui, po, sim, c_id, tl
--    from MDQueue
--    order by num_terms_left, heads
--    limit 1;
--
--    if heads = Null then
--      heads = '{}'::text[];
--    end if;
--
--    if tl = Null then
--      tl = '{}'::text[];
--    end if;
--
--
--    if found then
--  
--      delete from MDQueue where entry_id = id;
--  
--      if l > 0 then
--        term = tails[1];
--        tails = tails[2:l];
--  
--        if substring(term from 1 for 1) = '%'
--        then
--          insert into MDQueue (protein_id, expanded_terms, rest_of_terms, num_terms_left, 
--                               num_expjobs, num_iterations, num_steps_sim, user_id, policy,
--                               simulator, config_id, trj_list)
--          (select p_id, heads || symbol, tails, l-1, ne, ni, ns, ui, po, sim, c_id, tl
--           from AminoAcids 
--           where property similar to term
--           order by heads || symbol
--          );
--        else
--          insert into MDQueue (protein_id, expanded_terms, rest_of_terms, num_terms_left, 
--                               num_expjobs, num_iterations, num_steps_sim, user_id, policy,
--                               simulator, config_id, trj_list)
--          (select p_id, heads || symbol, tails, l-1, ne, ni, ns, ui, po, sim, c_id, tl
--           from AminoAcids 
--           where symbol similar to term
--           order by heads || symbol
--          );
--        end if;
--      end if;
--    end if;
--
--  end loop;

--  return p_id || ',' || array_to_string(heads, '') || ',' || ne || ',' || 
--         ni ||',' || ns || ',' || ui || ',' || po || ',' || sim || ',' || c_id;
--end;
--$$ language plpgsql;


drop function if exists controlparams_upsert(text, text);
create or replace function controlparams_upsert(n text, v text)
returns void as $$
begin
  update ControlParams set param_value = v where param_name = n;
  if not found then
    insert into ControlParams values (n,v);
  end if;
end;
$$ language plpgsql;




drop function if exists controlparams_getjson();
create or replace function controlparams_getjson()
returns text as $$
declare
  lines text[];
  res   text;
begin
  select array_agg('"' || param_name || '": "' || param_value || '"') 
  into lines 
  from ControlParams;

  res = E'{\n  ' || array_to_string(lines, E',\n  ') || E'\n}';

  return res;
end;
$$ language plpgsql;



drop function if exists subspaces_insert(bigint, bigint, text, int, int, text);
create or replace function subspaces_insert(c_id bigint, p_id bigint, gh text, gp int, u_id int, po text)
returns integer as $$
declare
  res int;
begin
  lock table Subspaces in exclusive mode;
  insert into Subspaces (config_id, protein_id, gmdhost, gmdport, t_start, user_id, policy) 
  values (c_id, p_id, gh ,gp, now(), u_id, po);
  select max(subspace_id) into res from Subspaces;
  return res;
end;
$$ language plpgsql;

drop function if exists subspaces_activate(s_id int);
create or replace function subspaces_activate(s_id int)
returns void as $$
declare
  e_id     int;
  iter     int;
begin
  lock table Subspaces in exclusive mode;
  lock table ExpJobs in exclusive mode;
  lock table RL_History in exclusive mode;

  insert into RL_History (subspace_id, state_id, action_id) values (s_id, null, null);
  select coalesce(max(iter_id),0) into iter
  from RL_History
  where subspace_id = s_id;

  for e_id in (select expjob_id from ExpJobs where subspace_id = s_id) loop
    perform expjobs_activate(e_id);
  end loop;
  perform schedule_straight_md_jobs(s_id);
end;
$$ language plpgsql;

drop function if exists subspaces_deactivate(s_id int);
create or replace function subspaces_deactivate(s_id int)  
returns void as $$
declare
  e_id int;
begin
  lock table Subspaces in exclusive mode;
  lock table ExpJobs in exclusive mode;

  for e_id in (select expjob_id from ExpJobs where subspace_id = s_id) loop
    perform expjobs_deactivate(e_id);
  end loop;
end;
$$ language plpgsql;


drop function if exists expjobs_insert(bigint, int[]);
create or replace function expjobs_insert(s_id bigint, w_desc int[])
returns integer as $$
declare
  res int;
  d_id bigint;
begin

  select expdesc_getId(w_desc) into d_id;
  lock table ExpJobs in exclusive mode;
  insert into ExpJobs (subspace_id, desc_id) values (s_id, d_id);
  select max(expjob_id) into res from ExpJobs;
  return res;
end;
$$ language plpgsql;

drop function if exists expjobs_insert(bigint);
create or replace function expjobs_insert(s_id bigint)
returns integer as $$
declare
  res int;
begin

  lock table ExpJobs in exclusive mode;
  insert into ExpJobs (subspace_id, desc_id) values (s_id, Null);
  select max(expjob_id) into res from ExpJobs;
  return res;
end;
$$ language plpgsql;



drop function if exists expjobs_insert(bigint, int, float, int);
create or replace function expjobs_insert(s_id bigint,  
                                          ns_sim int, ts_pico float,
                                          tsav_freq float,
                                          n int,
                                          simu text)
returns integer as $$
declare
  res int;
  d_id bigint;
  w_desc int[];
  num_features int;

begin
  select count(*) into num_features
  from ProteinFeatures P, Subspaces S
  where S.subspace_id = s_id and
        P.config_id   = S.config_id and
        P.protein_id  = S.protein_id;


  select array_agg(B.b order by B.b) into w_desc
  from generate_series(1,num_features) as B;


  select expdesc_getId(w_desc) into d_id;
  lock table ExpJobs in exclusive mode;

  insert into ExpJobs (subspace_id, desc_id, nsteps_sim, timestep_pico, trjsave_freq, 
                       num_iterations, total_num_iterations, simulator) 
  values (s_id, d_id, ns_sim, ts_pico, tsav_freq, n, n, simu);

  select max(expjob_id) into res from ExpJobs;
  return res;
end;
$$ language plpgsql;


create or replace function expjobs_insert(s_id bigint)
returns integer as $$
declare
  res int;
  d_id bigint;
  w_desc int[];
  num_features int;
begin
  select count(*) into num_features
  from ProteinFeatures P, Subspaces S
  where S.subspace_id = s_id and
        P.config_id   = S.config_id and 
        P.protein_id  = S.protein_id;


  select array_agg(B.b order by B.b) into w_desc
  from generate_series(1,num_features) as B;

  select expjobs_insert(s_id, w_desc) into res;
  return res;
end;
$$ language plpgsql;





drop function if exists expdescs_getId(w_desc int[]);
create or replace function expdesc_getId(w_desc int[])
returns integer as $$
declare
  res int;
  d   int[];
begin
  select array_agg(B.b order by B.b) into d from unnest(w_desc) as B;

  lock table ExpDescs in exclusive mode;
  select desc_id into res from ExpDescs where exp_desc = d;
  if not found then
    insert into ExpDescs (exp_desc) values (d);
    select max(desc_id) into res from ExpDescs;
  end if;

  return res;
end;
$$ language plpgsql;

drop function if exists expjobs_activate(e_id int);
create or replace function expjobs_activate(e_id int)
returns void as $$
declare
begin
  if exists (select 1 from ExpJobs where expjob_id = e_id and active) then
    return;
  end if;

  update ExpJobs set active = true where expjob_id = e_id;

end;
$$ language plpgsql;


drop function if exists expjobs_deactivate(e_id int);
create or replace function expjobs_deactivate(e_id int)
returns void as $$
declare
begin
  if exists (select 1 from ExpJobs where expjob_id = e_id and not active) then
    return;
  end if;

  update ExpJobs set active = false where expjob_id = e_id;                                                                         
end;
$$ language plpgsql;










