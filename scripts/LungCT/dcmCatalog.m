function T = dcmCatalog(varargin)
% Catalogs desired DICOM header values
% Inputs (optional):
%       extout = extension of files to save
%       pathout = directory for saving converted files

% Input: filt = string or cell array of strings for file name filters
%               e.g. '*.dcm' or {'*.1','*.dcm'}

opath = pwd;
filt = {'*.1','*.dcm','','*.IMA'};
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

dcmtags = {'PatientName','StudyID','StudyDate','PatientID','SeriesNumber',...
    'AccessionNumber','AcquisitionNumber','StudyDescription','SeriesDescription',...
    'ScanOptions','FilterType','ConvolutionKernel','ManufacturerModelName','Manufacturer','KVP',...
    'XRayTubeCurrent','ExposureTime','Exposure','SliceThickness',...
    'Rows','Columns','PixelSpacing'};
ntags = length(dcmtags);
tvars = [{'Directory'},dcmtags(1:(end-1)),{'dy','dx','Slices'}];
vtypes = {'cellstr',...
    'cellstr','cellstr','cellstr','cellstr','uint32',...
    'uint32','uint16','cellstr','cellstr',...
    'cellstr','cellstr','cellstr','cellstr','cellstr','double',...6
    'double','double','double','double',...
    'uint16','uint16',...
    'double','double','uint16'};
nvars = length(tvars);

% Find folders containing DICOMs
fprintf('Finding DICOM folders ...\n');
[D,F] = dirtree(opath,filt); % D: cell array of strings (directories)
% Remove DICOMDIR
ind = cellfun(@(x)(length(x)==1)&&strcmp(x{1},'DICOMDIR'),F);
F(ind) = [];
D(ind) = [];
ndir = length(D);

T = table('Size',[ndir,nvars],'VariableTypes',vtypes,'VariableNames',tvars);

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
        info = dicominfo(tname,'UseDictionaryVR',true);
        for ifld = 1:ntags
            ttag = dcmtags{ifld};
            if isfield(info,ttag)

                if strcmp(ttag,'PatientName')
                    T.PatientName{idir} = info.PatientName.FamilyName;
                elseif strcmp(ttag,'PixelSpacing')
                    T.dx(idir) = info.PixelSpacing(1);
                    T.dy(idir) = info.PixelSpacing(2);
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
        T.Directory{idir} = D{idir};
        T.Slices(idir) = length(F{idir});
    end
end
T(~Tflag,:) = []; % Remove rows for non-DICOM
fprintf('\n');

% Choose where to save results:
% [fname,opath] = uiputfile('*.csv','Save DICOM Results','DICOMcatalog.csv');

if ischar(opath)
    % Save results as CSV:
    writetable(T,fname);
end