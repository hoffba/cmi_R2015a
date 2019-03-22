function job = batch_fluxLungReg(username,catalog_name)
% Input:
%   catalog_name = name of text file cataloging lung data to process



% Set Elastix options
% elxbindir = sprintf('/scratch/turbo/umms-cgalban/Elastix');
elxfname = fullfile('/nfs/turbo/umms-cgalban/Elastix',{'ElastixParameters_1.txt','ElastixParameters_2.txt','ElastixParameters_3.txt'});
opts = {'f_clim',[-inf,0],...
        'm_clim',[-inf,0],...
        'f_dilateN',5*ones(1,3),...
        'f_filtN',ones(1,3),...
        'm_filtN',ones(1,3),...
        'jac',true};
    
% Read catalog file:    
T = readtable(catalog_name,'Format','%u%s%s');
uT = unique(T(:,1:2));
nuT = size(uT,1);

% Loop over cases:
for i = 1:nuT
    
    % Find registration files:
    fnames = pathFromFlux(table2array(T( (T.Var1==uT.Var1(i)) & strcmp(T.Var2,uT.Var2(i)) ,3)));
    
    % Determine results directories:
    [basepath,basename] = fileparts(fnames{1});
    % Desired output directory:
    odir = fullfile(basepath,[basename,'_elxreg']);
    % Temporarily use scratch storage for speed:
    workingdir = fullfile('/scratch/cgalban_fluxod',username,[basename,'_elxreg']);
    
    f_name = fnames{~cellfun(@isempty,strfind(fnames,'_exp.mhd'))};
    fv_name = fnames{~cellfun(@isempty,strfind(fnames,'_exp_label.mhd'))};
    m_name = fnames{~cellfun(@isempty,strfind(fnames,'_ins.mhd'))};
    mv_name = fnames{~cellfun(@isempty,strfind(fnames,'_ins_label.mhd'))};
    
    fprintf('Starting batch: %s\n',basename);
    
    job(i) = batch(@autoRegLung,0,[{odir,workingdir,elxfname,f_name,fv_name,m_name,mv_name},opts]);
    job(i).Tag = basename;
end

% Monitor jobs until all are finished:
running = true(nuT,1);
while any(running)
    pause(300);
    job_state = {job.State}';
    job_finished = strcmp(job_state,'finished');
    job_failed = strcmp(job_state,'failed');
    
    str = [{job.Tag};job_state'];
    fprintf('\nJob status:\n');
    fprintf('%s : %s\n',str{:});
    
    running = ~(job_finished|job_failed);
end



