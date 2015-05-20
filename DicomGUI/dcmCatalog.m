function C = dcmCatalog(extout,pathout)
% Catalogs desired DICOM header values
% Inputs (optional):
%       extout = extension of files to save
%       pathout = directory for saving converted files

% Input: filt = string or cell array of strings for file name filters
%               e.g. '*.dcm' or {'*.1','*.dcm'}

% if nargin==0
    filt = {'*.1','*.dcm',''};
% elseif ischar(filt)
%     filt = {filt};
% end
%possible parfor loop for increased processing speed?

hdrstr = {'PatientID','StudyID','StudyDate','SeriesNumber','AcquisitionNumber',...
    'StudyDescription','SeriesDescription',...
    'ScanOptions','FilterType','ConvolutionKernel','KVP',...
    'SliceThickness','Rows','Columns','PixelSpacing'};
nfld = length(hdrstr);
C = [{'Directory'},hdrstr(1:(end-1)),{'dy','dx','Slices'}];

% Find folders containing DICOMs
opath = uigetdir(pwd,'Select folder containing all DICOMS:');
if opath~=0
    hw = waitbar(0,'Finding DICOM folders ...');
    D = dirtree(opath,filt); % D: cell array of strings (directories)
    ndir = length(D);

    % Loop over all DICOM directories
    for idir = 1:ndir
        waitbar(idir/ndir,hw,['Processing folder ',num2str(idir),...
                         ' of ',num2str(ndir)]);
        fnames = [];
        for ifilt = 1:length(filt)
            tfnames = dir(fullfile(D{idir},filt{ifilt}));
            if isempty(filt{ifilt})
                tfnames(~cellfun(@isempty,strfind({tfnames(:).name},'.'))) = [];
                tfnames([tfnames(:).isdir]) = [];
            end
            fnames = cat(1,fnames,tfnames);
        end
        nf = length(fnames);
        if nf~=0
            fields = {}; clear('dinfo');
            for ifn = nf:-1:1
                tinfo = dicominfo(fullfile(D{idir},fnames(ifn).name));
                tfields = fieldnames(tinfo);
                newfields = setdiff(tfields,fields); enf = ~isempty(newfields);
                misfields = setdiff(fields,tfields); emf = ~isempty(misfields);
                if (ifn==nf)
                    fields = tfields;
                elseif (enf || emf)
                    % Need to adjust info structure fields to match
                    if emf
                        for k = 1:length(misfields)
                            tinfo.(misfields{k}) = [];
                        end
                    end
                    if enf
                        for k = 1:length(newfields)
                            dinfo(ifn).(newfields{k}) = [];
                        end
                        fields = [fields;newfields]; %#ok<*AGROW>
                        dinfo = orderfields(dinfo,tinfo);
                    else
                        tinfo = orderfields(tinfo,dinfo);
                    end
                end
                dinfo(ifn) = tinfo;
            end
            % Find unique acquisition numbers
            if isfield(dinfo,'AcquisitionNumber') ...
                    && ~any(cellfun(@isempty,{dinfo(:).AcquisitionNumber}))
                acqn = [dinfo(:).AcquisitionNumber];
                acqu = unique(acqn);
            else
                acqn = 0; acqu = 0;
            end
            % Store relevant parameters for each acquisition
            for j = 1:length(acqu)
                ii = find(acqn==acqu(j));
                tC = cell(1,nfld);
                for ifld = 1:nfld
                    if isfield(dinfo,hdrstr{ifld})
                        tC{ifld} = dinfo(ii(1)).(hdrstr{ifld});
                    end
                end
                if isempty(tC{end})
                    dd = nan(1,2);
                else
                    dd = tC{end};
                end
                C = [C;D(idir),tC(1:(end-1)),{dd(1),dd(2),length(ii)}];
            end
        end
        clear info
    end
    delete(hw);
    [fname,path] = uiputfile('*.csv','Save Catalog Data:');
    if fname~=0
        cmi_csvwrite(fullfile(path,fname),C);
    end
end
