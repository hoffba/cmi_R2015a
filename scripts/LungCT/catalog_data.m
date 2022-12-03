function T = catalog_data(path)

T = [];
vars = {'CaseNumber',             'uint16',       false;...
        'Tag',                    'cellstr',      false;...
        'UMlabel',                'cellstr',      false;...
        'PatientName',            'cellstr',      true;...
        'StudyDate',              'cellstr',      true;...
        'SeriesNumber',           'uint32',       true;...
        'ConvolutionKernel',      'cellstr',      true;...
        'SliceThickness',         'double',       true;...
        'Slices',                 'uint16',       false;...
        'DataPath',               'cellstr',      false;...
        'DataType',               'cellstr',      false;...
        'StudyID',                'cellstr',      true;...
        'PatientID',              'cellstr',      true;...
        'AccessionNumber',        'uint32',       true;...
        'AcquisitionNumber',      'uint16',       true;...
        'StudyDescription',       'cellstr',      true;...
        'SeriesDescription',      'cellstr',      true;...
        'ScanOptions',            'cellstr',      true;...
        'FilterTypes',            'cellstr',      true;...
        'ManufacturerModelName',  'cellstr',      true;...
        'Manufacturer',           'cellstr',      true;...
        'KVP',                    'double',       true;...
        'XRayTubeCurrent',        'double',       true;...
        'ExposureTime',           'double',       true;...
        'Exposure',               'double',       true;...
        'Rows',                   'uint16',       true;...
        'Columns',                'uint16',       true;...
        'dy',                     'double',       false;... % Actually looks for PixelSpacing
        'dx',                     'double',       false};
nvar = size(vars,1);
Tdef = table('Size',[1,nvar],'VariableNames',vars(:,1)','VariableTypes',vars(:,2)');

% Search for image files
filtstr = {'mhd','nii','nii.gz'};
for i = 1:numel(filtstr)
    fn = dir(fullfile(path,'**',['*.',filtstr{i}]));
    for j = 1:numel(fn)
        % Extract UMlabel and Date
        tok = regexp(fn.name,'(.*)_(\d{8})[_.]','tokens');
        if ~isempty(tok)
            umlabel = tok{1};
            datstr = tok{2};
        else
            umlabel = extractBefore(fn.name,'.');
            datstr = 'Unknown';
        end
        tT = Tdef;
        tT.UMlabel = umlabel;
        tT.PatientName = umlabel;
        tT.Studydate = datstr;
        tT.DataPath = fullfile(fn.folder,fn.name);
        tT.DataType = filtstr{i};
        T = [T;tT];
    end
end

% Search for DICOMs
filtstr = {'1.*','*.1','*.dcm','','*.IMA'};
[D,F] = dirtree(opath,filtstr);
for i = 1:numel(D)
    fname = fullfile(D{i},F{i}{1});
    if isdicom(fname)
        info = dicominfo(fname);
        tT = Tdef;
        tT.DataPath = D{i};
        tT.DataType = 'DICOM';
        tT.Slices = numel(F{i});
        if isfield(info,'PixelSpacing')
            tT.dx = info.PixelSpacing(1);
            tT.dy = info.PixelSpacing(2);
        end
        for j = 1:nvar
            if isfield(info,vars{i,1})
                tT.(vars{i,1}) = info.(vars{i,1});
            end
        end
        
        % Determine UM label:
        if isempty(info.PatientName) || strcmp(tT.UMlabel,'Anonymous')
            tT.UMlabel = info.PatientID;
        else
            tT.UMlabel = info.PatientName;
        end
        tT.UMlabel = regexprep(tT.UMlabel,' ','_');
        
        T = [T;tT];
    end
end

% Sort data, and determine groupings
T = sortrows(T,{'StudyID','PatientName','StudyDate','SeriesNumber'});
% Find case groupings
[~,~,ugroups_ic] = unique(strcat(T.StudyID,T.PatientName,T.StudyDate));
T.CaseNumber = ugroups_ic;

writetable(T,fullfile(path,'Pipeline_catalog.csv'));
