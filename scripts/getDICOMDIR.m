function T = getDICOMDIR(dcm_path)

fn = dir(dcm_path);
fn(1:2) = [];
N = numel(fn);

vtypes = {'cellstr','cellstr','cellstr','double','cell'};
vnames = {'StudyDate','StudyID','PatientName','SeriesNumber','DicomList'};
T = table('Size',[0,numel(vnames)],'VariableTypes',vtypes,'VariableNames',vnames);

fprintf('%5.1f%%',0);
for i = 1:numel(fn)
    if ~mod(i,round(N/1000))
        fprintf('\b\b\b\b\b\b%5.1f%%',i/N*100);
    end
    if isdicom(fullfile(dcm_path,fn(i).name))
        info = dicominfo(fullfile(dcm_path,fn(i).name));
        ind = find(strcmp(T.StudyDate,info.StudyDate) & strcmp(T.StudyID,info.StudyID) ...
            & strcmp(T.PatientName,info.PatientName.FamilyName) & info.SeriesNumber==info.SeriesNumber,1);
        [~,tname] = fileparts(info.Filename);
        if isempty(ind)
            % Add new series to catalog
            T = [T;{{info.StudyDate},{info.StudyID},{info.PatientName.FamilyName},info.SeriesNumber,{{tname}}}];
        else
            % Add DICOM file to existing series list
            T.DicomList{ind}{end+1,1} = tname;
        end
    end
end
                                                                                                                
% save(fullfile(dcm_path,'DICOMDIR.mat'),'T');

