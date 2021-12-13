function Nnodes = GL_sbatch(fname,jobname,GLstr,username,fn_inputs,inputs,nc,per_mem,per_time,N)
% Script for determining GL system parameters and saving SBATCH input file
% Inputs:
%       fname:      full file name for saving .sh file
%       jobname:    name of the job
%       GLstr:      name of the GL function being run
%       username:   UM uniquename of user
%       fn_inputs:  full file name of function inputs
%       nf:         number of cases to be run
%       pmem:       conservative estimate for required memory per case
%       ptime:      conservative estimate for required processing time per case
%       N:          (optional) number of requested nodes to array

% Run on Standard Node:
part_str = 'standard';
mx_mem = 180;  % 180GB max memory for standard node
mx_cores = 35; % 36 cores max on standard node, leave one for coordinating batch jobs
mx_time = 60*24*10; % Limit the runtime to 10 days

% Determine number of nodes requeseted
if nargin==7 && ~isempty(N) && round(N)>0
    Nnodes = round(N);
    N = true;
else
    Nnodes = 1;
    N = false;
end

% Limit on memory and number of cores:
cores = min([ floor(mx_mem / per_mem) , nc/Nnodes , mx_cores ]);

% Determine limitations due to walltime:
N_time = per_time * ceil( nc / ( cores * Nnodes ));
if N_time > mx_time
    if N
        warning(['Requested number of nodes results in a potential walltime of %u hrs,\n',...
            ' which is longer than the accepted limits (%u hrs) on Great Lakes.'],N_time/60,mx_time/60);
    else
        Nnodes = ceil( N_time / mx_time );
    end
end

% Calculate final requirements:
nc = nc / Nnodes;
cores = min([ floor(mx_mem / per_mem) , nc , mx_cores ]);
mem = per_mem * cores;
walltime = per_time * ceil(nc/cores);

% Print system parameters:
fprintf('Nodes : %u\nCores : %u\nMemory : %.1f GB\nWalltime : %s\n',...
    Nnodes,cores,mem,duration(0,walltime,0));

% Save Matlab inputs:
inputs.Nnodes = Nnodes;
save(fn_inputs,'-struct','inputs');

% Write SBATCH file
opts = {'partition',part_str, 'cores',cores, 'mem',mem, 'walltime',walltime,...
    'username',username, 'inputs_fname',fn_inputs};
if Nnodes > 1
    opts = [opts,{'array',Nnodes}];
end
write_sbatch(fname,jobname,GLstr,opts{:});
