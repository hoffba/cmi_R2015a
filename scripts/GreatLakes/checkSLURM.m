function [stat,str] = checkSLURM(ncases,part,mem,nnodes)

stat = true;
str = '';

% Determine what type of partition to use
partition_stats = {'standard',160;...
                   'largemem',1500};
if strcmp(part,'auto')
    if mem > 20
        opts.Partition = 'largemem';
    else
        opts.Partition = 'standard';
    end
end
opts.MaxMem = partition_stats{strcmp(opts.Partition,partition_stats(:,1)),2};

    % Calculate number of nodes
    % Max cores are adjusted smalled than max to make sure jobs can find an available partition
    switch part 
        case 'standard'
            core_max = 30;
        case 'largemem'
            core_max = 24;
        otherwise
            core_max = 20;
    end 
    cores = min(min(opts.Niter,floor(opts.MaxMem/opts.ProcessMemory))+1,core_max);
    walltime = opts.ProcessTime*ceil(opts.Niter/(cores-1));
    if isnan(opts.Nodes)
        opts.Nodes = ceil(walltime/(14*24*60)); % 14 days max walltime
    end

    % Recalculate taking into account #nodes
    opts.Niter = ceil(opts.Niter/opts.Nodes);
    cores = min(min(opts.Niter,floor(opts.MaxMem/opts.ProcessMemory))+1,36);
    mem = opts.ProcessMemory * (cores-1);
    walltime = opts.ProcessTime*ceil(opts.Niter/(cores-1));
    if walltime > (14*24*60) % 14 day max walltime
        warning('Processes may not complete.\nEstimated time = %s\n',...
            datestr(duration(0,walltime,0),'dd-HH:MM:SS'));
    end
