function p = getSLURM

C = {'jobid',       'SLURM_JOBID';...
     'jobname',     'SLURM_JOB_NAME';...
     'submitdir',   'SLURM_SUBMIT_DIR';...
     'submithost',  'SLURM_SUBMIT_HOST';...
     'nodelist',    'SLURM_JOB_NODELIST';...
     'arraytaskid', 'SLURM_ARRAY_TASK_ID';...
     'jobpartition','SLURM_JOB_PARTITION';...
     'nnodes',      'SLURM_JOB_NUM_NODES';...
     'ntasks',      'SLURM_NTASKS';...
     'procpernode', 'SLURM_TASKS_PER_NODE';
     'taskspernode','SLURM_NTASKS_PER_NODE';...
     'cpupertask',  'SLURM_CPUS_PER_TASK';
     'priority',    'SLURM_PRIO_PROCESS';
     'user',        'SLURM_JOB_USER';
     'hostname',    'SLURM_SUBMIT_HOST'};

p = struct([]);
for i = 1:size(C,1)
    val = getenv(C{i,2});
    valn = str2double(val);
    if isnan(valn)
        p.(C{i,1}) = val;
    else
        p.(C{i,1}) = valn;
    end
end