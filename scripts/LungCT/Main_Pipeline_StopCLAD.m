function j = Main_Pipeline_StopCLAD(dcmpath,svpath)

%% Find files
if exist(fullfile(dcmpath,'DICOMcatalog.csv'),'file')
    cases = catalog_select_2(dcmpath);
else
    C = dcmCatalog(dcmpath);
    cases = catalog_select_2(C);
end

%% Loop over cases
ncases = numel(cases);
for i = 1:ncases
    
    basename = sprintf('%s_%s',cases(i).PatientName,cases(i).StudyDate);
    
    procdir = fullfile(svpath,basename);
    
    
    %% Start pipeline batch job
    fprintf('%s - starting Main_Pipeline_sub as batch job: #',basename);
    j(i) = batch(@Main_Pipeline_sub,1,[{cases(i).Scans.Directory},{procdir}]);
    fprintf('%u\n',j(i).ID);
end
