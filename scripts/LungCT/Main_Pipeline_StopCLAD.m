function j = Main_Pipeline_StopCLAD(dcmpath,svpath)

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
%     fprintf('%s - starting Main_Pipeline_sub as batch job: #',basename);
%     j(i) = batch(@Main_Pipeline_sub,1,[{basename},{cases(i).Scans.Directory},{procdir}]);
%     fprintf('%u\n',j(i).ID);

    Main_Pipeline_sub(basename,cases(i).Scans.Directory,procdir);
end
