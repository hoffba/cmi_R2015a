function j = Main_Pipeline_select(varargin)
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

cluster_profile = {};
flag = false;
opts = struct('dcmpath','',...
              'sv_path','',...
              'quickreg',false);

if nargin
    % Find directory inputs
    pathi = find(cellfun(@(x)ischar(x) && isfolder(x),varargin),2);
    if numel(pathi)
        opts.dcmpath = varargin{pathi(1)};
    end
    if numel(pathi)>1
        opts.sv_path = varargin{pathi(2)};
    end
    if ismember('quickreg',varargin)
        opts.quickreg = true;
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
                fprintf('Could not find cluster profile: %s\n',mls_fname);
                cluster_profile = {};
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
    
    if flag
        % Start pipeline in command window
        fprintf('%s - starting Main_Pipeline_sub case #%u\n',basename,i);
        Main_Pipeline_sub(basename,cases(i).Scans.Directory,procdir,opts);
    else
        % Start pipeline batch job
        fprintf('%s - starting Main_Pipeline_sub as batch job: #',basename);
        j(i) = batch(@Main_Pipeline_sub,1,[{basename},{cases(i).Scans.Directory},{procdir},{opts}],cluster_profile{:});
        fprintf('%u\n',j(i).ID);
    end
end
