function j = Main_Pipeline_select(dcmpath,svpath,varargin)
% Start CT lung pipeline script (Main_Pipeline_sub) with GUI selection for
% cases based on DCMcatalog.csv
%
% Main_Pipeline_select( dcmpath, svpath )
%       * runs cases using full reg and default local cluster settings
% Main_Pipeline_select( dcmpath, svpath, 'quickreg' )
%       * skips last step of registration for quick results
% Main_Pipeline_select( dcmpath, svpath, 'serial' )
%       * runs cases serially in the command line (for debugging)
% Main_Pipeline_select( dcmpath, svpath, 'serialbatch' )
%       * runs cases serially in batch (1-worker)
% Main_Pipeline_select( dcmpath, svpath, 'quickreg', 'serialbatch' )

quickreg = false;
cluster_profile = {};
flag = false;

if ismember('quickreg',varargin)
    quickreg = true;
end
if ismember('serial',varargin)
    flag = true;
elseif ismember('serialbatch',varargin)
    cluster_profile = {'Profile','cmi_1'};
    if ~ismember('cmi_1',parallel.clusterProfiles)
        mls_fname = fullfile(fileparts(which('cmi')),'cmi_1.mlsettings');
        if exist(mls_fname,'file')
            fprintf('Importing cluster profile: %s\n',mls_fname);
            parallel.importProfile(mls_fname);
        else
            cluster_profile = {};
        end
    end
end

%% Find files
if ischar(dcmpath)
    if isfolder(dcmpath)
        dcmpath = fullfile(dcmpath,'DICOMcatalog.csv');
    end
    cases = catalog_select_2(dcmpath);
else
    error('Invalid input.')
end

%% Loop over cases
ncases = numel(cases);
for i = 1:ncases
    
    basename = sprintf('%s_%s',cases(i).UMlabel,cases(i).StudyDate);
    procdir = fullfile(svpath,basename);
    
    if flag
        % Start pipeline in command window
        fprintf('%s - starting Main_Pipeline_sub case #%u\n',basename,i);
        Main_Pipeline_sub(basename,cases(i).Scans.Directory,procdir,quickreg);
    else
        % Start pipeline batch job
        fprintf('%s - starting Main_Pipeline_sub as batch job: #',basename);
        j(i) = batch(@Main_Pipeline_sub,1,[{basename},{cases(i).Scans.Directory},{procdir},{quickreg}],cluster_profile{:});
        fprintf('%u\n',j(i).ID);
    end
end
