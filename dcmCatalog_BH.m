function C = dcmCatalog_BH
% Catalogs desired DICOM header values
% Inputs (optional):
%       extout = extension of files to save
%       pathout = directory for saving converted files

% Input: filt = string or cell array of strings for file name filters
%               e.g. '*.dcm' or {'*.1','*.dcm'}

filt = {'*.1','*.dcm',''};

dcmtags = {'PatientID','StudyID','StudyDate','SeriesNumber','AcquisitionNumber',...
    'StudyDescription','SeriesDescription','ScanOptions',...
    'SliceThickness','Rows','Columns','PixelSpacing'};
nfld = length(dcmtags);
C = [{'Directory'},dcmtags(1:(end-1)),{'dy','dx','Slices'}];

% Find folders containing DICOMs
opath = uigetdir(pwd,'Select folder containing all DICOMS:');
if opath~=0
    hw = waitbar(0,'Finding DICOM folders ...');
    [D,F] = dirtree(opath,filt); % D: cell array of strings (directories)
    % Remove DICOMDIR
    ind = cellfun(@(x)(length(x)==1)&&strcmp(x{1},'DICOMDIR'),F);
    F(ind) = [];
    D(ind) = [];
    ndir = length(D);

    % Loop over all DICOM directories
    C = [{'Directory'},dcmtags(1:(end-1)),{'dy','dx','Slices'};cell(ndir,nfld+3)];
    for idir = 1:ndir
        waitbar(idir/ndir,hw,['Processing folder ',num2str(idir),...
                         ' of ',num2str(ndir)]);
        info = dicominfo_BH(fullfile(D{idir},F{idir}{1}),'UseDictionaryVR',true);
        tC = cell(1,nfld);
        for ifld = 1:nfld
            if isfield(info,dcmtags{ifld})
                tC{ifld} = info.(dcmtags{ifld});
            end
        end
        if isempty(tC{end})
            dd = nan(1,2);
        else
            dd = tC{end};
        end
        C(idir+1,:) = [D(idir),tC(1:(end-1)),{dd(1),dd(2),length(F{idir})}];
    end
    delete(hw);
    
    % Choose where to save results:
    [fname,fpath] = uiputfile('*.csv','Save DICOM Results','DICOMcatalog.csv');
    
    if ischar(fpath)
        % Save results as CSV:
        cmi_csvwrite(fullfile(fpath,fname),C);
    end
end
