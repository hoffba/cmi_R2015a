function j = Main_Pipeline_select(varargin)
% Start CT lung pipeline script (Main_Pipeline_sub) with GUI selection for
% cases based on DCMcatalog.csv
%
% Main_Pipeline_select( dcm_path, save_path )
%       * runs cases using full reg and default local cluster settings
% Main_Pipeline_select( dcm_path, save_path, 'quickreg' )
%       * skips last step of registration for quick results
% Main_Pipeline_select( dcm_path, save_path, 'serial' )
%       * runs cases serially in the command line (for debugging)
% Main_Pipeline_select( dcm_path, save_path, 'serialbatch' )
%       * runs cases serially in batch (1-worker)
% Main_Pipeline_select( dcm_path, save_path, 'quickreg', 'serialbatch' )

cluster_profile = {};

% Processing flag: True means run serially on command window
flag = false;

opts = struct('dcm_path','',...
              'save_path','',...
              'quickreg',false);

%% Determine input options
if nargin
    opts.dcm_path = varargin{1};
    if nargin>1
        opts.save_path = varargin{2};
        if nargin>2
            opts.quickreg = ismember('quickreg',varargin);
            if ismember('serial',varargin)
                flag = true;
%             elseif ismember('serialbatch',varargin)
%                 flag = 2;
%                 cluster_profile = {'Profile','cmi_1'};
%                 if ~ismember('cmi_1',parallel.clusterProfiles)
%                     mls_fname = fullfile(fileparts(which('cmi')),'cmi_1.mlsettings');
%                     if exist(mls_fname,'file')
%                         fprintf('Importing cluster profile: %s\n',mls_fname);
%                         parallel.importProfile(mls_fname);
%                     else
%                         fprintf('Could not find cluster profile: %s\n',mls_fname);
%                         cluster_profile = {};
%                     end
%                 end
            end
        end
    end
end

%% Find files
[cases,opts] = catalog_select_3('opts',opts);

%% Loop over cases
ncases = numel(cases);
for i = 1:ncases
    
    basename = sprintf('%s_%s',cases(i).UMlabel,cases(i).StudyDate);
    procdir = fullfile(opts.save_path,basename);
    if flag % serial
        % Start pipeline in command window
        fprintf('%s - starting Main_Pipeline_sub case #%u\n',basename,i);
        Main_Pipeline_sub(basename,cases(i).Scans.Directory,procdir,opts);
    else
        if i==1
            myCluster = parcluster;
            myCluster.NumWorkers = opts.par_size;
        end
        fprintf('%s - starting Main_Pipeline_sub case #%d\n',basename,i);
        j(i) = batch(myCluster,@Main_Pipeline_sub,1,[{basename},{cases(i).Scans.Directory},{procdir},{opts}]);
    end
%         case 2 % serialbatch
%             % Start pipeline batch job
%             fprintf('%s - starting Main_Pipeline_sub as batch job: #',basename);
%             j(i) = batch(@Main_Pipeline_sub,1,[{basename},{cases(i).Scans.Directory},{procdir},{opts}],cluster_profile{:});
%             fprintf('%u\n',j(i).ID);
%     end
end
