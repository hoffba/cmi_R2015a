function T = catalog_data(searchpath,gflag)

if nargin<2
    gflag = false;
end

T = [];

if gflag
    hw = waitbar(0,'Searching for files ...');
end

% Search for DICOMs
filtstr = {'1.*','*.1','*.dcm','*.IMA'};
[D,F] = dirtree(searchpath,filtstr);
% Files to exclude (endsWith):
exlstr = {'.bmp'};

if isempty(D)
    % Search for image files
    filtstr = {'mhd','nii','nii.gz'};
    for i = 1:numel(filtstr)
        fn = dir(fullfile(searchpath,'**',['*.',filtstr{i}]));
        for j = 1:numel(fn)
            % Extract UMlabel and Date
            %   UMlabel becomes the first string before _
            %   StudyDate becomes the first 8-digit number after that

            tok = regexp(fn(j).name,'([^\.]+)(.*)','tokens');
            namestr = tok{1}{1};
            ext = tok{1}{2};

            t = getDICOMvars;
            t.SeriesDescription = {fn(j).name};
            nametok = strsplit(namestr,'_');
            ind = 1;
            t.PatientName = nametok(1);
            if numel(nametok)>1 && any(cellfun(@(x)~isempty(regexp(nametok{2},x,'once')),{'\d{8}','t\d+'}))
                t.StudyDate = nametok(2);
                ind = 2;
            end
            t.UMlabel = {strjoin(nametok(1:ind),'_')};
            nametok(1:ind) = [];
            if ~isempty(nametok)
                t.StudyDescription = {strjoin(nametok,'_')};
            end
            t.DataPath = {fullfile(fn(j).folder,fn(j).name)};
            t.DataType = filtstr(i);

            % Try and find tags
            tset = false;
            if contains(ext,'.exp.label.')
                t.Tag{1} = 'ExpLabel'; tset = true;
            elseif contains(ext,'.exp.')
                t.Tag{1} = 'Exp'; tset = true;
            elseif contains(ext,'.ins.label.')
                t.Tag{1} = 'InsLabel'; tset = true;
            elseif contains(ext,'.ins.')
                t.Tag{1} = 'Ins'; tset = true;
            end
            if ~tset
                TF = contains({'label','voi','seg'},nametok,'IgnoreCase',true);
                if any(strcmpi(nametok,'exp'))
                    if any(TF)
                        t.Tag = {'ExpLabel'};
                    else
                        t.Tag = {'Exp'};
                    end
                elseif any(strcmpi(nametok,'ins'))
                    if any(TF)
                        t.Tag = {'InsLabel'};
                    else
                        t.Tag = {'Ins'};
                    end
                end
            end

            T = [T;t];
        end
    end
else
    ND = numel(D);
    for i = 1:ND
        if gflag && isvalid(hw)
            waitbar(i/ND,hw,sprintf('Cataloging; %d / %d',i,ND));
        end
        % Exclude files with extension:
        for j = 1:numel(exlstr)
            F{i}(endsWith(F{i},exlstr{j})) = [];
        end

        dcm_flag = false;
        if ~isempty(F{i})
            fname = fullfile(D{i},F{i}{1});
            try
                dcm_flag = isdicom(fname);
            catch err
                disp(fname)
            end
        end
        if dcm_flag
            info = dicominfo(fname,'UseDictionaryVR',true);
            if all(isfield(info,{'StudyDate','SeriesNumber','AccessionNumber','AcquisitionNumber'}))
                [t,vname] = getDICOMvars(info.Modality);
                t.DataPath = D(i);
                t.DataType = {'DICOM'};
                t.Slices = numel(F{i});
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
end

if ~isempty(T)
    % Sort data, and determine groupings
    T = sortrows(T,{'StudyID','PatientName','StudyDate','SeriesNumber'});
    % Find case groupings
    [~,~,ugroups_ic] = unique(strcat(T.StudyID,T.PatientName,T.StudyDate));
    T.CaseNumber = ugroups_ic;

    writetable(T,fullfile(searchpath,'Data_Catalog.csv'));
end

if gflag && isvalid(hw)
    delete(hw);
end

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
        if ischar(t.(newvar{i})) || iscellstr(t.(newvar{i}))
            T.(newvar{i}) = repmat({''},Nr,1);
        elseif isstring(t.(newvar{i}))
            T.(newvar{i}) = repmat("",Nr,1);
        else
            T.(newvar{i}) = cast(nan(Nr,1),class(t.(newvar{i})));
        end
    end

