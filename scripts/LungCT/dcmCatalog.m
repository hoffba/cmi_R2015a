function T = dcmCatalog(varargin)
% Catalogs desired DICOM header values
% Inputs (optional):
%       extout = extension of files to save
%       pathout = directory for saving converted files

% Input: filt = string or cell array of strings for file name filters
%               e.g. '*.dcm' or {'*.1','*.dcm'}

opath = pwd;
filt = {'1.*','*.1','*.dcm','','*.IMA'};
fname = '';
if nargin
    for i = 1:nargin
        val = varargin{i};
        if ischar(val) && isfolder(val) % Folder selection
            opath = val;
        elseif ischar(val) && strcmp(val(end-3:end),'.csv')
            fname = val;
        elseif iscellstr(val) && any(contains(val,'*'))
            filt = val;
        end
    end
else
    % Use GUI to select folder:
    opath = uigetdir(pwd,'Select folder containing all DICOMS:');
    if ~opath
        return;
    end
end
if isempty(fname)
    fname = fullfile(opath,'DICOMcatalog.csv');
end
%         VariableName              VariableType    DICOM field flag
vnames = {'CaseNumber',             'uint16',       false;...
          'Tag',                    'cellstr',      false;...
          'UMlabel',                'cellstr',      false;...
          'PatientName',            'cellstr',      true;...
          'StudyDate',              'cellstr',      true;...
          'SeriesNumber',           'uint32',       true;...
          'ConvolutionKernel',      'cellstr',      true;...
          'SliceThickness',         'double',       true;...
          'Slices',                 'uint16',       false;...
          'Directory',              'cellstr',      false;...
          'StudyID',                'cellstr',      true;...
          'PatientID',              'cellstr',      true;...
          'AccessionNumber',        'uint32',       true;...
          'AcquisitionNumber',      'uint16',       true;...
          'StudyDescription',       'cellstr',      true;...
          'SeriesDescription',      'cellstr',      true;...
          'ScanOptions',            'cellstr',      true;...
          'FilterTypes',            'cellstr',      true;...
          'ManufacturerModelName',  'cellstr',      true;...
          'Monufacturer',           'cellstr',      true;...
          'KVP',                    'double',       true;...
          'XRayTubeCurrent',        'double',       true;...
          'ExposureTime',           'double',       true;...
          'Exposure',               'double',       true;...
          'Rows',                   'uint16',       true;...
          'Columns',                'uint16',       true;...
          'dy',                     'double',       true;... % Actually looks for PixelSpacing
          'dx',                     'double',       false};
nvars = size(vnames,1);

% Find folders containing DICOMs
fprintf('Finding DICOM folders ...\n');
[D,F] = dirtree(opath,filt); % D: cell array of strings (directories)
% Remove DICOMDIR
ind = cellfun(@(x)(length(x)==1)&&strcmp(x{1},'DICOMDIR'),F);
F(ind) = [];
D(ind) = [];
ndir = length(D);

T = table('Size',[ndir,nvars],'VariableTypes',vnames(:,2)','VariableNames',vnames(:,1)');

% Loop over all DICOM directories
fprintf('Processing folder (of %u): ',ndir);
Tflag = true(ndir,1);
nchar = 0;
for idir = 1:ndir
    ntxt = sprintf('%u',idir);
    fprintf([repmat('\b',1,nchar),'%s'],ntxt);
    nchar = numel(ntxt);
    
    tname = fullfile(D{idir},F{idir}{1});
    Tflag(idir) = isdicom(tname);
    if Tflag(idir)
        
        % Find DICOM fields of interest:
        info = dicominfo(tname,'UseDictionaryVR',true);
        for ifld = 1:nvars
            if vnames{ifld,3}
                ttag = vnames{ifld,1};
                if strcmp(ttag,'dy') && isfield(info,'PixelSpacing')
                    T.dx(idir) = info.PixelSpacing(1);
                    T.dy(idir) = info.PixelSpacing(2);
                elseif isfield(info,ttag)
                    if strcmp(ttag,'PatientName')
                        T.PatientName{idir} = info.PatientName.FamilyName;
                    else
                        val = info.(ttag);
                        if ischar(val)
                            if iscellstr(T.(ttag))
                                T.(ttag){idir} = val;
                            else
                                T.(ttag)(idir) = str2double(val);
                            end
                        elseif ~isempty(val)
                            T.(ttag)(idir) = val;
                        end
                    end

                end
            end
        end
        
        % Set other info based on DICOM info
        T.Directory{idir} = D{idir};
        T.Slices(idir) = length(F{idir});
        
    end
end
T(~Tflag,:) = []; % Remove rows for non-DICOM
fprintf('\n');

% Determine a label for our use
T.UMlabel = T.PatientName;
% Try setting empty/anonymouse names using PatientID instead
ind = cellfun(@(x)isempty(x)||strcmp(x,'Anonymous'),T.UMlabel);
T.UMlabel(ind) = T.PatientID(ind);

% Sort and group scans:
% Sort by identifiers
T = sortrows(T,{'StudyID','PatientName','StudyDate','SeriesNumber'});
% Find case groupings
[~,~,ugroups_ic] = unique(strcat(T.StudyID,T.PatientName,T.StudyDate));
T.CaseNumber = ugroups_ic;
% For cases with no identifiers, name by group number
ind = cellfun(@isempty,T.UMlabel);
T(ind).UMlabel = cellfun(@(x)sprintf('%05.0f',x),ugroups_ic(ind),'UniformOutput',false);



% Choose where to save results:
% [fname,opath] = uiputfile('*.csv','Save DICOM Results','DICOMcatalog.csv');

if ischar(opath)
    % Save results as CSV:
    writetable(T,fname);
end