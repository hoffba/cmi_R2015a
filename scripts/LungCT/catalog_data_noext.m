function T = catalog_data_noext(path)

T = [];

% Search for DICOMs
fprintf('Finding DICOM subfolders...\n')
D = dir('**\*');
D = unique({D([D.isdir]).folder});

N = numel(D);
ndig = floor(log10(N))+1;
fmt_str = sprintf('%%%dd/%d',ndig,N);
fprintf(['Cataloging DICOM series... ',fmt_str],0);
fmt_str = [repmat('\b',1,2*ndig+1),fmt_str];
for i = 1:numel(D)
    fprintf(fmt_str,i);
    F = dir(D{i});
    fname = '';
    for fi = 1:numel(F)
        try
            tfname = fullfile(F(fi).folder,F(fi).name);
            if isdicom(tfname)
                fname = tfname;
                break
            end
        end
    end
    if ~isempty(fname)
        info = dicominfo(fname,'UseDictionaryVR',true);
        if all(isfield(info,{'StudyDate','SeriesNumber','AccessionNumber','AcquisitionNumber'}))
            [t,vname] = getDICOMvars(info.Modality);
            t.DataPath = D(i);
            t.DataType = {'DICOM'};

            % Estimate number of slices
            [~,fn] = fileparts(fname);
            str = extractBefore(fn,digitsPattern);
            t.Slices = nnz(startsWith({F.name},str));

            if isfield(info,'PixelSpacing')
                t.dx = info.PixelSpacing(1);
                t.dy = info.PixelSpacing(2);
            end
            for j = 1:numel(vname)
                if isfield(info,vname{j})
                    if strcmp(vname{j},'PatientName')
                        val = info.PatientName.FamilyName;
                    else
                        val = info.(vname{j});
                    end
                    if ischar(val)
                        val = {val};
                    end
                    if ~isempty(val)
                        t.(vname{j}) = val;
                    end
                end
            end

            % Determine UM label:
            if isfield(info,'PatientName') && ~isempty(info.PatientName) && ~strcmp(info.PatientName.FamilyName,'Anonymous')
                t.UMlabel = {info.PatientName.FamilyName};
            elseif isfield(info,'PatientID')
                t.UMlabel = {info.PatientID};
            elseif isfield(info,'StudyID')
                t.UMlabel = info.StudyID;
            else
                t.UMlabel = sprintf('SUBJ%5.0f',i);
            end
            t.UMlabel = regexprep(t.UMlabel,' ','_');

            T = addTableRow(T,t);
        end
    end
end

% Sort data, and determine groupings
T = sortrows(T,{'StudyID','PatientName','StudyDate','SeriesNumber'});
% Find case groupings
[~,~,ugroups_ic] = unique(strcat(T.StudyID,T.PatientName,T.StudyDate));
T.CaseNumber = ugroups_ic;

catname = fullfile(path,'Data_Catalog.csv');
fprintf('\nSaving catalog to: %s\n',catname);
writetable(T,catname);

function [t,vnames] = getDICOMvars(modstr)
%       Catalog value             Value class
vnames = {'CaseNumber',             'uint16';...
          'Tag',                    'cellstr';...
          'UMlabel',                'cellstr';...
          'PatientName',            'cellstr';...
          'StudyDate',              'cellstr';...
          'SeriesNumber',           'uint32';...
          'SliceThickness',         'double';...
          'Slices',                 'uint16';...
          'DataPath',               'cellstr';...
          'DataType',               'cellstr';...
          'StudyID',                'cellstr';...
          'PatientID',              'cellstr';...
          'AccessionNumber',        'cellstr';...
          'AcquisitionNumber',      'uint16';...
          'StudyDescription',       'cellstr';...
          'SeriesDescription',      'cellstr';...
          'ScanOptions',            'cellstr';...
          'Manufacturer',           'cellstr';...
          'ManufacturerModelName',  'cellstr';...
          'Rows',                   'uint16';...
          'Columns',                'uint16';...
          'dy',                     'double';... % Actually looks for PixelSpacing
          'dx',                     'double'};
if nargin
    switch modstr
        case 'CT'
            vnames = [vnames;...
                    {'ConvolutionKernel',  'cellstr';...
                     'FilterTypes',        'cellstr';...
                     'KVP',                'double';...
                     'XRayTubeCurrent',    'double';...
                     'ExposureTime',       'double';...
                     'Exposure',           'double'}];
        case 'MR'
            vnames = [vnames;...
                    {'ScanningSequence',   'cellstr';...
                     'MRAcquisitionType',  'cellstr';...
                     'RepetitionTime',     'double';...
                     'EchoTime',           'double'}];
    end
end
t = table('Size',[1,size(vnames,1)],'VariableTypes',vnames(:,2)','VariableNames',vnames(:,1)');
vnames = vnames(:,1)';

function T = addTableRow(T,t)
if isempty(T)
    T = t;
else
    % Find any new variable to add to table
    T = addVars(T,t);
    % Add missing variables to input to match the output
    t = addVars(t,T);
    % Add new
    idx = zeros(1,size(T,2));
    for i = 1:numel(idx)
        idx(i) = find(strcmp(T.Properties.VariableNames{i},t.Properties.VariableNames),1);
    end
    T = [T;t(:,idx)];
end

function T = addVars(T,t)
Nr = size(T,1);
T_vnames = T.Properties.VariableNames;
t_vnames = t.Properties.VariableNames;
newvar = t_vnames(~ismember(t_vnames,T_vnames));
for i = 1:numel(newvar)
    if ischar(t.(newvar{i}))
        T.(newvar{i}) = repmat({''},Nr,1);
    else
        T.(newvar{i}) = cast(nan(Nr,1),class(t.newvar{i}));
    end
end
