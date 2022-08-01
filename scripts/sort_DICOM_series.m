function T = sort_DICOM_series(dcm_path)

fn_dcmdir = fullfile(dcm_path,'DICOMDIR.mat');
% if exist(fn_dcmdir,'file')
%     T = load(fn_dcmdir);
%     T = T.T;
% else
    fn = dir(dcm_path);
    fn(1:2) = [];
    N = numel(fn);

    vtypes = {'cellstr','cellstr','cellstr','double','cell'};
    vnames = {'StudyDate','StudyID','PatientName','SeriesNumber','DicomList'};
    T = table('Size',[0,numel(vnames)],'VariableTypes',vtypes,'VariableNames',vnames);

    fprintf('Reading DICOM info: %4.1f%%',0);
    for i = 1:numel(fn)
        if ~mod(i,round(N/1000))
            fprintf('\b\b\b\b\b%4.1f%%',i/N*100);
        end
        if isdicom(fullfile(dcm_path,fn(i).name))
            info = dicominfo(fullfile(dcm_path,fn(i).name));
            ind = find(strcmp(T.StudyDate,info.StudyDate) & strcmp(T.StudyID,info.StudyID) ...
                & strcmp(T.PatientName,info.PatientName.FamilyName) & T.SeriesNumber==info.SeriesNumber,1);
            dirname = sprintf('%s%s_%s_%u',info.StudyID,info.PatientName.FamilyName,info.StudyDate,info.SeriesNumber);
            tname = flip(extractBefore(flip(info.Filename),filesep));
            if isempty(ind)
                % Add new series to catalog
                T = [T;{{info.StudyDate},{info.StudyID},{info.PatientName.FamilyName},info.SeriesNumber,{{tname}}}];
            else
                % Add DICOM file to existing series list
                T.DicomList{ind}{end+1,1} = tname;
            end
            if ~isfolder(dirname)
                mkdir(dirname);
            end
            if ~exist(fullfile(dirname,tname),'file')
                copyfile(fullfile(dcm_path,fn(i).name),fullfile(dirname,tname));
            end
        end
    end
    fprintf('\n');

    save(fn_dcmdir,'T');

% end