--
-- Molecular dynamics DB schema
--
-- Author: Yanif Ahmad, 04/01/2012
--

drop view if exists BenchmarkView;
drop view if exists AtomView;
drop view if exists SS_Histograms;
drop view if exists ProteinNumDims;
drop view if exists ExpJobsView;
drop view if exists JobDurationsView;

drop table if exists ConfigVDWNonBonded;
drop table if exists ConfigDihedrals;
drop table if exists ConfigAngles;
drop table if exists ConfigBonds;
drop table if exists ConfigAtoms;

drop table if exists Restraints;


drop table if exists ProteinFeatures;
drop table if exists ProteinDihedrals;
drop table if exists ProteinAngles;
drop table if exists ProteinBonds;
drop table if exists ProteinAtoms;
drop table if exists ProteinSizes;

drop table if exists AtomPositions cascade;
drop table if exists AtomVelocities;
drop table if exists Energies;
drop table if exists ConformationSpace;
drop table if exists ConformationReservoir;
drop table if exists ReservoirStats;
drop table if exists FinishedTrajectories;
drop table if exists CurrIterConformationSpace;

drop table if exists HistogramRegionScores;
drop table if exists Histograms;
drop table if exists BiasedHistograms;
drop table if exists Buckets;
drop table if exists HistogramSpec;
drop table if exists Partitions;
drop table if exists PartitionedHistograms;
drop table if exists IterHistograms;
drop table if exists IterHistogramsHI;

drop table if exists PartitionedHistogramsCache;
drop table if exists IterHistogramsHICache;
drop table if exists IterHistogramsCache;
drop table if exists CacheTimestamps;
drop table if exists CacheLogs;

drop table if exists TrajectoryPIDs;
drop table if exists MachinePerformance;
drop table if exists TrajectoryPerformance;
drop table if exists WorkerPerformance;
drop table if exists DatabasePerformance;
drop table if exists MachineInfo;

drop table if exists SystemStats;

drop table if exists ErrorMessages;


drop table if exists HistoBinAssignments;


drop table if exists Trajectories cascade;

drop table if exists Proteins cascade;
drop table if exists ExpJobs;
drop table if exists ExpDescs cascade;
drop table if exists Subspaces cascade;

drop table if exists AtomTypeIDs cascade;

drop table if exists SimConfigs cascade;


create table SimConfigs (
  config_id         serial,
  simulator         text, -- (Charmm/Amber/...)
  ff                text, -- (Amber94/Charmm97/...)
  primary key       (config_id)
);



create table Proteins (
  protein_id        serial,
  protein_cname     text, -- e.g. AA for Alanine Dipeptide
  protein_name      text, -- human readable description of the protein system e.g. 'alanine_dipeptide'
  protein_pdb_id    text,
  protein_fn        text,
  primary key       (protein_id)
);
create index Proteins_cname_index on Proteins (protein_cname);

create table AtomTypeIDs (
  atom_type_id serial,
  atom_type    text,
  primary key (atom_type_id)
);
create index AtomTypeIDs_atom_type_Index on AtomTypeIDs (atom_type);

create table Subspaces (
  subspace_id      serial,
  user_id          bigint,
  config_id        bigint references SimConfigs(config_id),
  protein_id       bigint references Proteins(protein_id),
  gmdhost          text,
  gmdport          int,
  blocking         boolean   default true,
  num_waits        int       default null,
  timeout_secs     int       default 0,
  t_start          timestamp with time zone,
  status           text,
  policy           text,
  primary key      (subspace_id)
);

create table ExpDescs (
  desc_id          serial,
  exp_desc         int[],
  primary key      (desc_id)
);
create index ExpDesc_cg_desc_Index on ExpDescs (exp_desc);

create table ExpJobs (
  expjob_id            serial,
  subspace_id          bigint references Subspaces(subspace_id),
  desc_id              bigint references ExpDescs(desc_id),
  active               boolean default false,
  blocked              boolean default false,
  nsteps_sim           int     default pow(10,5),
  timestep_pico        float   default 0.002,
  trjsave_freq         int     default pow(10,2),
  num_iterations       int,
  total_num_iterations int,
  simulator            text    default 'namd',
  primary key          (expjob_id)
);
create index ExpJobs_subspace_id_Index on ExpJobs (subspace_id);


---------------------------------------------------------------------------------------------
-- Config Specific
---------------------------------------------------------------------------------------------

create table ConfigAtoms (
  config_id    bigint references SimConfigs(config_id),
  atom_type_id bigint references AtomTypeIDs(atom_type_id),
  residue_name text,
  atom_name    text,
  charge       float,
  primary key (config_id, residue_name, atom_type_id, atom_name)
);

create table ConfigBonds (
  config_id         bigint references SimConfigs(config_id),
  atom_type_id1     bigint references AtomTypeIDs(atom_type_id),
  atom_type_id2     bigint references AtomTypeIDs(atom_type_id),
  bond_const        float,
  bond_length       float,
  primary key  (config_id, atom_type_id1, atom_type_id2)
);

create table ConfigAngles (
  config_id         bigint references SimConfigs(config_id),
  atom_type_id1     bigint references AtomTypeIDs(atom_type_id),
  atom_type_id2     bigint references AtomTypeIDs(atom_type_id),
  atom_type_id3     bigint references AtomTypeIDs(atom_type_id),
  angle_const       float,
  angle             float,
  primary key  (config_id, atom_type_id1, atom_type_id2, atom_type_id3)
);

create table ConfigDihedrals (
  config_id         bigint references SimConfigs(config_id),
  atom_type_id1     bigint references AtomTypeIDs(atom_type_id),
  atom_type_id2     bigint references AtomTypeIDs(atom_type_id), 
  atom_type_id3     bigint references AtomTypeIDs(atom_type_id),
  atom_type_id4     bigint references AtomTypeIDs(atom_type_id),
  improper          boolean,
  force_const       float, -- Kchi in degrees
  n                 float, -- multiplicity
  delta             float, -- phase shift
  primary key    (config_id, atom_type_id1, atom_type_id2, atom_type_id3, atom_type_id4)
);

create table ConfigVDWNonBonded (
  config_id         bigint references SimConfigs(config_id),
  atom_type_id1     bigint references AtomTypeIDs(atom_type_id),
  atom_type_id2     bigint references AtomTypeIDs(atom_type_id),
  vdw_rmin          double precision,
  vdw_eps           double precision,
  primary key (config_id, atom_type_id1, atom_type_id2)
);

---------------------------------------------------------------------------------------------
-- Protein specific topological information which is also config specific
-- because different force fields use different coarse graining models and that
-- affects the topology of the protein.
---------------------------------------------------------------------------------------------

create table ProteinAtoms (
  config_id         bigint references SimConfigs(config_id),
  protein_id        bigint references Proteins(protein_id),
  atom_id           int,
  atom_type         text,
  atom_type_id      bigint references AtomTypeIDs(atom_type_id),
  atom_name         text,
  residue_id        int,
  residue_name      text,
  segment_name      text,
  primary key (config_id, protein_id, atom_id)
);

create table ProteinBonds (
  config_id         bigint references SimConfigs(config_id),
  protein_id        bigint references Proteins(protein_id),
  atom_id1          int, --references ProteinAtoms(atom_id),
  atom_id2          int, --references ProteinAtoms(atom_id),
  primary key  (config_id, protein_id, atom_id1, atom_id2)
);

create table ProteinAngles (
  config_id         bigint references SimConfigs(config_id),
  protein_id        bigint references Proteins(protein_id),
  atom_id1          int, --references ProteinAtoms(atom_id),
  atom_id2          int, --references ProteinAtoms(atom_id),
  atom_id3          int, --references ProteinAtoms(atom_id),
  primary key  (config_id, protein_id, atom_id1, atom_id2, atom_id3)
);

create table ProteinDihedrals (
  config_id         bigint references SimConfigs(config_id),
  protein_id        bigint references Proteins(protein_id),
  atom_id1          int, --references ProteinAtoms(atom_id),
  atom_id2          int, --references ProteinAtoms(atom_id),
  atom_id3          int, --references ProteinAtoms(atom_id),
  atom_id4          int, --references ProteinAtoms(atom_id),
  improper          boolean,
  primary key  (config_id, protein_id, atom_id1, atom_id2, atom_id3, atom_id4)
);

create table ProteinSizes (
  config_id         bigint references SimConfigs(config_id),
  protein_id        bigint references Proteins(protein_id),
  protein_size      int,
  primary key (config_id, protein_id)
);




create table ProteinFeatures (
  config_id         bigint references SimConfigs(config_id),
  protein_id        bigint references Proteins(protein_id),
  atom_id1          int, --references ProteinAtoms(atom_id),
  atom_id2          int, --references ProteinAtoms(atom_id),
  atom_id3          int, --references ProteinAtoms(atom_id),
  atom_id4          int, --references ProteinAtoms(atom_id),
  feature_id        bigint,
  primary key  (config_id, protein_id, atom_id1, atom_id2, atom_id3, atom_id4)
);
create index ProteinFeatures_idIndex on ProteinFeatures (feature_id);

create or replace view ProteinNumDims as
select config_id,protein_id,max(feature_id) as num_dims 
from proteinfeatures group by config_id,protein_id;


drop table if exists RL_Actions cascade;
create table RL_Actions (
  action_id     serial,
  peak_ratio    float,
  primary key   (action_id)
);

drop table if exists RL_History;
create table RL_History (
  iter_id     serial,
  subspace_id bigint references Subspaces(subspace_id),
  state_id    bigint,
  action_id   bigint references RL_Actions(action_id),
  acc_score   float,
  bumpiness   float,
  primary key (iter_id)
);

---------------------------------------------------------------------------------------------
-- MSM Related
---------------------------------------------------------------------------------------------

drop table if exists MSM_EnergyWellsXYZ;
create table MSM_EnergyWellsXYZ (
  subspace_id  bigint,
  well_id      bigint,
  ref_id       int,
  x            float,
  y            float,
  z            float,
  primary key  (subspace_id, well_id)
);

drop table if exists MSM_EnergyWellsPhiPsi;
create table MSM_EnergyWellsPhiPsi (
  subspace_id  bigint,
  well_id      bigint,
  ref_id       int,
  feature_id   int,
  centroid     float,
  span         float, -- min = centroid - span/2; max = centroid + span/2
  primary key  (subspace_id, well_id, ref_id, feature_id)
);


drop table if exists MSM_Transitions cascade;
create table MSM_Transitions (
  subspace_id  bigint,
  expjob_id    bigint,
  tran_id      serial,
  well_id1     bigint,
  well_id2     bigint,
  trj_id1      bigint,
  trj_id2      bigint
);
create index MSM_Transitions_index on MSM_Transitions (subspace_id, well_id1);

drop view if exists MSM;
create or replace view MSM as
select subspace_id, well_id1, well_id2, count(*) 
from msm_transitions 
where well_id2 is not null 
group by subspace_id, well_id1, well_id2 
order by subspace_id, well_id1, well_id2;


---------------------------------------------------------------------------------------------
-- RL_History
---------------------------------------------------------------------------------------------

drop table if exists RL_History;
create table RL_History (
  iter_id     serial,
  subspace_id int,
  state_id    int,
  action_id   int,
  acc_score   float,
  bumpiness   float,
  primary key (iter_id)
);


---------------------------------------------------------------------------------------------
-- Trajectory dependent and incrementally maintained
---------------------------------------------------------------------------------------------

create table Trajectories (
  trj_id         serial, -- (Also used as job_id)
  expjob_id      bigint references ExpJobs(expjob_id),
  biased         boolean default false,
  biased_feature bigint,
  iter_id        bigint references RL_History(iter_id),
  prev_trj_id    bigint,
  tran_id        bigint,
  pre_existing   boolean default false,
  step_size_fs   int,
  primary key    (trj_id)
);
create index Traj_iter_index on Trajectories (iter_id);
create index Traj_expjob_index on Trajectories (expjob_id, trj_id);

create table TrajectoryPIDs (
  trj_id         bigint references Trajectories(trj_id),
  process_id     int,
  primary key    (trj_id)
);

create table TrajectoryPerformance (
  trj_id         bigint references Trajectories(trj_id),
  time           timestamp with time zone,
  cpu            float,
  mem            float,
  dr_kb_s        float,
  dw_kb_s        float,
  io             float,
  primary key    (trj_id, time)
);

create table DatabasePerformance (
  time           timestamp with time zone,
  cpu            float,
  mem            float,
  dr_kb_s        float,
  dw_kb_s        float,
  io             float,
  primary key    (time)
);

create table MachinePerformance (
  name           text,
  time           timestamp with time zone,
  cpu            float,
  mem            float,
  dr_kb_s        float,
  dw_kb_s        float,
  io             float,
  primary key    (name, time)
);


create table FinishedTrajectories (
  trj_id        bigint references Trajectories(trj_id),
  primary key   (trj_id)
);


create table AtomVelocities (
  trj_id        bigint references Trajectories(trj_id),
  t             int,
  atom_id       int, -- references ProteinAtoms(atom_id),
  x             float,
  y             float,
  z             float
);
create index AtomVelocities_trj_id_index on AtomVelocities (trj_id);



create table AtomPositions (
  trj_id        bigint references Trajectories(trj_id),
  t             int,
  atom_id       int, -- references ProteinAtoms(atom_id),
  x             float,
  y             float,
  z             float
);
create index AtomPositions_trj_id_index on AtomPositions (trj_id);

create table Energies (
  trj_id        bigint references Trajectories(trj_id),
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
  tempavg       float,
  primary key   (trj_id, t)
);

create table ConformationSpace (
  --point_id       in
  trj_id         bigint references Trajectories(trj_id),
  t              int,
  feature_id     bigint, --references ProteinFeatures(feature_id),
  phi_psi        float8,
  primary key    (trj_id, t, feature_id)
); 
create index ConformationSpace_trj_idIndex on ConformationSpace (trj_id);

create table CurrIterConformationSpace (
  trj_id         bigint references Trajectories(trj_id),
  t              int,
  feature_id     bigint, -- references ProteinFeatures(feature_id),
  phi_psi        float8,
  primary key    (trj_id, t, feature_id)
);


create table Restraints (
  trj_id         bigint references Trajectories(trj_id),
  force          float,
  min            float,
  period         float
);
---------------------------------------------------------------------------------------------
-- Statistics for conformation space
---------------------------------------------------------------------------------------------
-- Sarana said comment this out
--create table Buckets (
--  num_dims       int,
--  feature_id     bigint references ProteinFeatures(feature_id),
--  bucket_id      int,
--  bucket_start   float8,
--  bucket_end     float8,
--  primary key    (num_dims, feature_id, bucket_start, bucket_end)
--);
--create index Buckets_id_index on Buckets (num_dims, feature_id, bucket_id);


create table HistogramSpec (
  spec_id             int,
  num_dims            int,
  buckets_per_dim     int[],
  regions_per_dim     int[],
  range_starts        float8[],
  range_ends          float8[]
);


create table BiasedHistograms (
  trj_id      bigint references Trajectories(trj_id),
  id_array    int[],
  region_id   int[],

  n           int default 0,
  s           int default 0
);
create index BiasedHistograms_index on BiasedHistograms (trj_id, id_array);


create table Histograms (
  trj_id      bigint references Trajectories(trj_id),
  id_array    int[],
  region_id   int[],
  
  n           int default 0,
  s           int default 0
);
create index Histograms_index on Histograms (trj_id, id_array);
create index Histogram_region_index on Histograms (trj_id,region_id);

create or replace view SS_Histograms as
select S.subspace_id,id_array,sum(n) as n,sum(s) as s
from Histograms H, trajectories T, ExpJobs E, SubSpaces S
where H.trj_id = T.trj_id and
      T.expjob_id = E.expjob_id and
      E.subspace_id = S.subspace_id
group by S.subspace_id,id_array;


create table HistogramRegionScores (
  trj_id      bigint references Trajectories(trj_id),
  region_id   int[],
  n           int,
  s           int,
  score       float
);
create index HistogramRegionScores_index on HistogramRegionScores (trj_id, region_id);




create table HistoBinAssignments (
  conf_id      bigint,
  iter_id      bigint references RL_History(iter_id),
  trj_id       bigint references Trajectories(trj_id),
  t            int,
  e            float8,
  bucket_id    int[],
  primary key  (conf_id)
);
create index HistoBinAssignments_index on HistoBinAssignments (iter_id,t);


create table PartitionedHistograms (
  iter_id      bigint references RL_History(iter_id),
  partition_id int,
  bucket_id    int[],
  n            int,
  s            int,
  conf_id      bigint references HistoBinAssignments(conf_id),
  e            float8,
  primary key  (iter_id, partition_id, bucket_id)
);

create table PartitionedHistogramsCache (
  iter_id      bigint references RL_History(iter_id),
  partition_id int,
  bucket_id    int[],
  n            int,
  s            int,
  conf_id      bigint references HistoBinAssignments(conf_id),
  e            float8,
  primary key  (iter_id, partition_id, bucket_id)
);

create table IterHistograms (
  iter_id     bigint references RL_History(iter_id),
  bucket_id   int[],
  n           int,
  s           int,
  conf_id     bigint references HistoBinAssignments(conf_id),
  e           float8,
  primary key (iter_id, bucket_id)
);

create table IterHistogramsHI (
  iter_id     bigint references RL_History(iter_id),
  bucket_id   int[],
  n           int,
  s           int,
  conf_id     bigint references HistoBinAssignments(conf_id),
  e           float8,
  primary key (iter_id, bucket_id)
);



create table IterHistogramsCache (
  iter_id     bigint references RL_History(iter_id),
  bucket_id   int[],
  n           int,
  s           int,
  conf_id     bigint references HistoBinAssignments(conf_id),
  e           float8,
  primary key (iter_id, bucket_id)
);

create table IterHistogramsHICache (
  iter_id     bigint references RL_History(iter_id),
  bucket_id   int[],
  n           int,
  s           int,
  conf_id     bigint references HistoBinAssignments(conf_id),
  e           float8,
  primary key (iter_id, bucket_id)
);


create table CacheTimestamps (                                                                                    
  iter_id     bigint references RL_History(iter_id),
  ts          timestamp with time zone,
  access_type text,
  note        text,
  primary key (iter_id)
);


create table CacheLogs (
  iter_id     bigint references RL_History(iter_id),
  ts          timestamp with time zone,
  access_type text,
  note        text,
  hit         boolean
);



---------------------------------------------------------------------------------------------
-- Sampling Control
---------------------------------------------------------------------------------------------

drop table if exists Benchmark;
drop table if exists GPUs;
drop table if exists Workers;
drop table if exists MDJobParams;
drop table if exists HistogramStats;
drop table if exists ControlParams;
drop table if exists AminoAcids;
drop table if exists MDQueue;
drop table if exists ClusterConfig;


create table ClusterConfig (
  gw_host      text,
  num_workers  int,
  task_id      text
);

create table Workers (
  gw_sn        serial,
  gw_host      text,
  gw_name      text,
  process_id   int,
  is_ongoing   boolean default true,
    -- Set it to fault if you want the worker to terminate after n_jobs_left jobs.
  n_jobs_left  int default 0,
    -- if the worker is not ongoing then n_jobs_left indicates
    -- the number of jobs it has to do before it can terminate.
  alive_signal timestamp with time zone,
  primary key  (gw_sn)
);

create table WorkerPerformance (
  gw_name        text,
  time           timestamp with time zone,
  cpu            float,
  mem            float,
  dr_kb_s        float,
  dw_kb_s        float,
  io             float,
  primary key    (gw_name, time)
);
--Sarana said we dont need it so ....
--create table MDJobParams (
--  trj_id       serial,
--  param_name   text,
--  param_value  text,
--  primary key  (trj_id, param_name)
--);

create table ControlParams (
  param_name   text,
  param_value  text,
  primary key  (param_name)
);


create table Benchmark (
  trj_id          bigint references Trajectories(trj_id),
  gw_name         text,
  gpu_worker      boolean,
  hostname        text,
  job_type        text,
  md_progress     int,
  ntime_steps     int,
  t1              timestamp with time zone,
  t2              timestamp with time zone,
  t3              timestamp with time zone,
  primary key     (trj_id)
);


create table HistogramStats (
  id                 serial,
  subspace_id        bigint references Subspaces(subspace_id),
  desc_id            bigint references ExpDescs(desc_id),
  curr_size          bigint,
  max_order_of_mag   float,
  primary key        (subspace_id,desc_id,id)
);

drop type if exists MDQEntry cascade;
create type MDQEntry as (
  protein_id         bigint,
  expanded_terms     text[],
  rest_of_terms      text[],
  num_terms_left     int,
  num_expjobs        int,
  num_iterations     int,
  num_steps_sim      int,
  trj_save_freq      int,
  timestep_pico      float,
  user_id            bigint,
  policy             text,
  simulator          text,
  config_id          int,
  trj_fn             text,
  psf_fn             text,
  well_fn            text
);


create table MDQueue (
  entry_id           serial,
  entry              MDQEntry,
  primary key        (entry_id)
);



create table AminoAcids (
  symbol   text,
  name     text,
  property text
);


create table GPUs (
  device_id   integer,
  gw_name     text,
  occupied    boolean
);

create table SystemStats (
  iter_id     bigint references RL_History(iter_id),
  n_frames    bigint,
  n_positions bigint,
  n_bins      bigint
);

create table ErrorMessages (
  ts          timestamp with time zone,
  place       text,
  message     text
);

create table MachineInfo (
  host_name   text,
  key         text,
  value       text
);

---------------------------------------------------------------------------------------------
-- Database Init
---------------------------------------------------------------------------------------------



insert into GPUs (device_id, occupied) 
select a.device_id, false as occupied 
from (select generate_series(0,7) as device_id)a;


--
--select d.a,d.b,variance(c) as var
--from (select s.a, r.b,sum(n) as c
--      from histogram h,
--      generate_series(1,20) as s(a),
--      generate_series(1,20) as r(b)
--      group by id_array[a],id_array[b],a,b)d
--group by d.a,d.b order by var;
--

-- select trj_id, array_agg(phi_psi order by feature_id) from ConformationSpace where t = 1 group by trj_id order by trj_id;

-- select S.trj_id, S.p as assigned, R.p as actual from StartingPoints S, (select trj_id, array_agg(phi_psi order by feature_id) p from ConformationSpace where t = 1 group by trj_id order by trj_id offset 2) R where S.trj_id = R.trj_id order by trj_id;


-- select c, count(*) as count_of_counts from (select distinct region_id[1:15],count(*) as c from histogram group by region_id[1:15])a group by c order by c;


create or replace view AtomView as
select atom_name,residue_name,segment_name,x,y,z 
from AtomPositions A, ProteinAtoms P, Trajectories T, ExpJobs E, Subspaces S 
where A.atom_id = P.atom_id and 
      A.trj_id = T.trj_id and 
      T.expjob_id = E.expjob_id and 
      E.subspace_id = S.subspace_id and 
      P.config_id = S.config_id and 
      P.protein_id = S.protein_id;


---------------------------------------------------------------------------------------------
-- RL
---------------------------------------------------------------------------------------------

-- --Sarana said comment this out
-- --drop table if exists RL_States;
-- --create table RL_States (
-- --  subspace_id bigint references Subspaces(subspace_id),
-- --  iter_id     bigint references RL_History(iter_id),
-- --  state_id    int,
-- --  region_id   int[],
-- --  n           int default 0
-- --);
-- create index RL_States_index on RL_States (subspace_id,iter_id,state_id,region_id);
 
drop table if exists RL_QTable;
create table RL_QTable (
  subspace_id bigint references Subspaces(subspace_id),
  state_id    bigint, -- references RL_History(state_id),
  action_id   bigint references RL_Actions(action_id),
  q_value     float,
  primary key (subspace_id,state_id,action_id)
);



---------------------------------------------------------------------------------------------
-- Views
---------------------------------------------------------------------------------------------



create or replace view BenchmarkView as
select B.gw_name, E.expjob_id, S.subspace_id, T.trj_id, 
       P.protein_cname as protein_seq, t1, (t3-t1) as duration
from Proteins P, Subspaces S, ExpJobs E, Trajectories T, Benchmark B 
where P.protein_id = S.protein_id and 
      S.subspace_id = E.subspace_id and
      E.expjob_id = T.expjob_id and
      B.trj_id = T.trj_id;

create or replace view ExpJobsView as
select S.subspace_id, S.t_start, P.protein_cname, E.expjob_id, 
       E.nsteps_sim, E.timestep_pico, E.num_iterations 
from subspaces S, proteins P, expjobs E, simconfigs C 
where S.subspace_id = E.subspace_id and 
      C.config_id = S.config_id and 
      S.protein_id = P.protein_id;

create or replace view JobDurationsView as
select S.subspace_id,T.expjob_id, H.iter_id,T.trj_id, B.t1 as t_start, B.t3-B.t1 as duration
from SubSpaces S, RL_History H, Trajectories T, Benchmark B
where S.subspace_id = H.subspace_id and
      T.iter_id     = H.iter_id and
      B.trj_id      = T.trj_id and
      not T.biased
order by S.subspace_id, H.iter_id, T.expjob_id, T.trj_id;



