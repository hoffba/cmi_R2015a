function T = dcmCatalog_BH_CJG(varargin)
% Catalogs desired DICOM metadata
% Syntax:
%   T = dcmCatalog_BH_CJG;
%   T = dcmCatalog_BH_CJG(filt);
%   T = dcmCatalog_BH_CJG(opath); 
%   T = dcmCatalog_BH_CJG(opath,svname);
%   T = dcmCatalog_BH_CJG(opath,svname,filt);
%
% Inputs (optional):
%       filt = cellstr with filename filters
%       opath = directory for searching for dicoms and saving results
%       svname = filename for output catalog .csv file

opath = '';
svname = '';
filt = {'*.1','*.dcm',''};
if nargin
    if iscellstr(varargin{1})
        filt = varargin{1};
    else
        opath = varargin{1};
    end
end
if nargin>1
    svname = varargin{2};
end
if nargin==3
    filt = varargin{3};
end
p = inputParser;
p.addRequired('filt',@iscellstr);
p.addRequired('opath',@(x) ischar(x) && isfolder(x) );
p.addRequired('svname',@(x) ischar(x) && isempty(regexp(x, '[/\*:?"<>|]', 'once')) );
p.parse(filt,opath,svname);
filt = p.Results.filt;
opath = p.Results.opath;
svname = p.Results.svname;

if isempty(opath)
    opath = uigetdir(pwd,'Select folder containing all DICOMS:');
end

dcmtags = {'PatientName','StudyID','StudyDate','PatientID','SeriesNumber','AccessionNumber',...
    'AcquisitionNumber','StudyDescription','SeriesDescription','ImageType',...
    'ScanOptions','FilterType','ConvolutionKernel','Manufacturer','KVP',...
    'XRayTubeCurrent','ExposureTime','Exposure','SliceThickness','Rows','Columns','PixelSpacing'};
vtype = {'cellstr','cellstr','cellstr','cellstr','cellstr','double','cellstr',...
    'double','cellstr','cellstr','cellstr','cellstr','cellstr','cellstr','cellstr','double',...
    'double','double','double','double','double','double','double','double','double'};
nfld = length(dcmtags);

% Find folders containing DICOMs
if opath~=0
    hw = waitbar(0,'Finding DICOM folders ...');
    [D,F] = dirtree(opath,filt); % D: cell array of strings (directories)
    % Remove DICOMDIR
    ind = cellfun(@(x)(length(x)==1)&&strcmp(x{1},'DICOMDIR'),F);
    F(ind) = [];
    D(ind) = [];
    ndir = length(D);

    % Loop over all DICOM directories
    vnames = [{'Directory'},dcmtags(1:(end-1)),{'dy','dx','Slices'}];
    T = table('Size',[ndir,numel(vnames)],'VariableTypes',vtype,'VariableNames',vnames);
    flagx = false(ndir,1);
    for idir = 1:ndir
        waitbar(idir/ndir,hw,['Processing folder ',num2str(idir),...
                         ' of ',num2str(ndir)]);
        
        fname = fullfile(D{idir},F{idir}{1});
        if isdicom(fname)
            info = dicominfo(fname,'UseDictionaryVR',true);
            for ifld = 1:nfld
                fieldname = dcmtags{ifld};
                if isfield(info,fieldname)
                    if strcmp(fieldname,'PatientName')
                        T.(fieldname){idir} = info.PatientName.FamilyName;
                    elseif strcmp(fieldname,'PixelSpacing')
                        T.dy(idir) = info.PixelSpacing(1);
                        T.dx(idir) = info.PixelSpacing(2);
                    else
                        tval = info.(fieldname);
                        if isnumeric(tval)
                            if isempty(tval)
                                T.(fieldname)(idir) = nan;
                            else
                                T.(fieldname)(idir) = tval;
                            end
                        else
                            T.(fieldname){idir} = tval;
                        end
                    end
                end
            end
            T.Slices(idir) = length(F{idir});
            T.Directory(idir) = D(idir);
        else
            flagx(idir) = true;
        end
    end
    delete(hw);
    T(flagx,:) = [];
    
    % Choose where to save results:
    if isempty(svname)
        [svname,opath] = uiputfile('*.csv','Save DICOM Results','DICOMcatalog.csv');
    end
    
    if svname
        % Save results as CSV:
        writetable(T,fullfile(opath,svname));
    end
end
