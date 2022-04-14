function j = Main_Pipeline_select(dcmpath,svpath,varargin)

quickreg = false;
cluster_profile = {};

if ismember('quickreg',varargin)
    quickreg = true;
end
if ismember('serial',varargin)
    cluster_profile = {'Profile','cmi_1'};
    if ~ismember('cmi_1',parallel.clusterProfiles)
        mls_fname = fullfile(fileparts(which('cmi')),'cmi_1.mlsettings');
        if exist(mls_fname,'file')
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
    
    %% Start pipeline batch job
    fprintf('%s - starting Main_Pipeline_sub as batch job: #',basename);
    j(i) = batch(@Main_Pipeline_sub,1,[{basename},{cases(i).Scans.Directory},{procdir},{quickreg}],cluster_profile(:));
    fprintf('%u\n',j(i).ID);
end
