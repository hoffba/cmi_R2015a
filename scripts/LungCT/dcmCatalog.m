function T = dcmCatalog(varargin)
% Catalogs desired DICOM header values
% Inputs (optional):
%       extout = extension of files to save
%       pathout = directory for saving converted files

% Input: filt = string or cell array of strings for file name filters
%               e.g. '*.dcm' or {'*.1','*.dcm'}

opath = pwd;
filt = {'*.1','*.dcm',''};
if nargin
    % Folder selection
    ind = find(isfolder(varargin),1);
    if ~isempty(ind)
        opath = varargin{ind};
    end
    % Filter selection
    ind = find(cellfun(@(x)any(contains(x,'*')),varargin),1);
    if ~isempty(ind)
        filt = varargin{ind};
    end
else
    % Use GUI to select folder:
    opath = uigetdir(pwd,'Select folder containing all DICOMS:');
    if ~opath
        return;
    end
end

dcmtags = {'PatientName','StudyID','StudyDate','PatientID','SeriesNumber',...
    'AccessionNumber','AcquisitionNumber','StudyDescription','SeriesDescription',...
    'ScanOptions','FilterType','ConvolutionKernel','Manufacturer','KVP',...
    'XRayTubeCurrent','ExposureTime','Exposure','SliceThickness',...
    'Rows','Columns','PixelSpacing'};
ntags = length(dcmtags);
tvars = [{'Directory'},dcmtags(1:(end-1)),{'dy','dx','Slices'}];
vtypes = {'char',...
    'char','char','char','char','uint32',...
    'uint32','uint16','char','char',...
    'char','char','char','char','double',...
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
nchar = 0;
for idir = 1:ndir
    ntxt = sprintf('%u',idir);
    nchar = numel(ntxt);
    fprintf([repmat()],
    waitbar(idir/ndir,hw,['Processing folder ',num2str(idir),...
        ' of ',num2str(ndir)]);
    
    info = dicominfo(fullfile(D{idir},F{idir}{1}),'UseDictionaryVR',true);
    tC = cell(1,ntags);
    for ifld = 1:ntags
        
        if isfield(info,dcmtags{ifld})
            
            if strcmp(dcmtags{ifld},'PatientName')
                tC{ifld} = info.PatientName.FamilyName;
            else
                tC{ifld} = info.(dcmtags{ifld});
            end
        end
    end
    if isempty(tC{end})
        dd = nan(1,2);
    else
        dd = tC{end};
    end
    T(idir+1,:) = [D(idir),tC(1:(end-1)),{dd(1),dd(2),length(F{idir})}];
end
delete(hw);

% Choose where to save results:
[fname,fpath] = uiputfile('*.csv','Save DICOM Results','DICOMcatalog.csv');

if ischar(fpath)
    % Save results as CSV:
    cmi_csvwrite(fullfile(fpath,fname),T);
end