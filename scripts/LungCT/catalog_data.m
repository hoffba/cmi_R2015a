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
        'AccessionNumber',        'cellstr',      true;...
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

% Search for DICOMs
filtstr = {'1.*','*.1','*.dcm','','*.IMA'};
[D,F] = dirtree(path,filtstr);

if isempty(D)
    % Search for image files
    filtstr = {'mhd','nii','nii.gz'};
    for i = 1:numel(filtstr)
        fn = dir(fullfile(path,'**',['*.',filtstr{i}]));
        for j = 1:numel(fn)
            % Extract UMlabel and Date
            %   UMlabel becomes the first string before _
            %   StudyDate becomes the first 8-digit number after that
            umlabel = ''; datstr = '';
            tok = regexp(fn(j).name,'(?:re_)?([^\.]*)_([^_\.]*)\.','tokens');
            if ~isempty(tok)
                umlabel = tok{1}{1};
                datstr = tok{1}{2};
            end
            tT = Tdef;
            tT.SeriesDescription = {fn(j).name};
            tT.UMlabel = {umlabel};
            tT.PatientName = {umlabel};
            tT.StudyDate = {datstr};
            tT.DataPath = {fullfile(fn(j).folder,fn(j).name)};
            tT.DataType = filtstr(i);
            T = [T;tT];
        end
    end
else
    for i = 1:numel(D)
        fname = fullfile(D{i},F{i}{1});
        if isdicom(fname)
            info = dicominfo(fname,'UseDictionaryVR',true);
            if all(isfield(info,{'StudyDate','SeriesNumber','AccessionNumber','AcquisitionNumber'}))
                tT = Tdef;
                tT.DataPath = D(i);
                tT.DataType = {'DICOM'};
                tT.Slices = numel(F{i});
                if isfield(info,'PixelSpacing')
                    tT.dx = info.PixelSpacing(1);
                    tT.dy = info.PixelSpacing(2);
                end
                for j = 1:nvar
                    vname = vars{j,1};
                    if isfield(info,vname)
                        if strcmp(vname,'PatientName')
                            val = info.PatientName.FamilyName;
                        else
                            val = info.(vars{j,1});
                        end
                        if ischar(val)
                            val = {val};
                        end
                        if ~isempty(val)
                            tT.(vars{j,1}) = val;
                        end
                    end
                end

                % Determine UM label:
                if isfield(info,'PatientName') && ~isempty(info.PatientName) && ~strcmp(info.PatientName.FamilyName,'Anonymous')
                    tT.UMlabel = {info.PatientName.FamilyName};
                elseif isfield(info,'PatientID')
                    tT.UMlabel = {info.PatientID};
                elseif isfield(info,'StudyID')
                    tT.UMlabel = info.StudyID;
                else
                    tT.UMlabel = sprintf('SUBJ%5.0f',i);
                end
                tT.UMlabel = regexprep(tT.UMlabel,' ','_');

                T = [T;tT];
            end
        end
    end
end

% Sort data, and determine groupings
T = sortrows(T,{'StudyID','PatientName','StudyDate','SeriesNumber'});
% Find case groupings
[~,~,ugroups_ic] = unique(strcat(T.StudyID,T.PatientName,T.StudyDate));
T.CaseNumber = ugroups_ic;

writetable(T,fullfile(path,'Pipeline_catalog.csv'));
