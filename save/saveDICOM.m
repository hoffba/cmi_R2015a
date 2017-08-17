function status = saveDICOM(fname,img,label,fov,info,verbchk)
status = false;
if (nargin>=4)
    if nargin<6
        verbchk = true;
    end
    
    [d(2),d(1),ns,nv] = size(img);
    % First check that the file name is correct
    [pathstr, fname, ~] = fileparts(fname);
    if ~isempty(label) && ~isnumeric(label)
        if ischar(label)
            label = {label};
        end
        if ~iscellstr(label) || (length(label)~=nv)
            label = strcat(label(1),num2cell(num2str((1:nv)')));
        end
    end
    
    % Auto-scale to int16:
    mxscale = 2^15-1;
    imax = squeeze(max(max(max(img,[],1),[],2),[],3));
    imin = squeeze(min(min(min(img,[],1),[],2),[],3));
    sthk = fov(3)/ns;
    spos = (sthk - fov(3))/2;
    warning('off','images:dicomwrite:inconsistentIODAndCreateModeOptions');
    
    % Determine DICOM header tags:
    if (nargin<5)
        info = {};
    end
    info = parseinfo(info);
    info.PixelSpacing = fov(1:2)./d;
    info.SpacingBetweenSlices = sthk;
    info.SliceThickness = sthk;
    
    for ivec = 1:nv
        info.AcquisitionNumber = ivec;
        info.RescaleIntercept = imin(ivec);
        info.RescaleSlope = (imax(ivec) - imin(ivec))/mxscale;
        if info.RescaleSlope == 0
            info.RescaleSlope = 1;
        end
        info.SeriesDescription = label{ivec};
        info.ImagePositionPatient = [0,0,spos];
        
        if nv>1
            fnout = [fname,'_',num2str(ivec)];
        else
            fnout = fname;
        end
        
        if ~exist(fullfile(pathstr,fnout),'dir')
            mkdir(fullfile(pathstr,fnout));
        end
        if verbchk, hw = waitbar(0,'Saving DICOM slices ...','Name','Saving DICOM image ...'); end
        for islc = 1:ns
            info.ImagePositionPatient(3) = info.ImagePositionPatient(3) + sthk;
            info.SliceLocation = info.ImagePositionPatient(3);
            info.InstanceNumber = islc;
            timg = int16((img(:,:,islc,ivec)-info.RescaleIntercept)/info.RescaleSlope);
            stat = dicomwrite(timg,fullfile(pathstr,fnout,sprintf([fnout,'_s%04u.dcm'],islc)),info);
            status = status && isempty(stat.BadAttribute) ...
                && isempty(stat.MissingCondition) ...
                && isempty(stat.MissingData) ...
                && isempty(stat.SuspectAttribute);
            if verbchk, waitbar(islc/ns,hw,['Slc ',num2str(info.ImagePositionPatient(3))]); end
        end
        if verbchk, close(hw); end
    end
end

function info = parseinfo(info)
info = [fieldnames(info)';struct2cell(info)'];
ind = strncmp('Private_',info(1,:),8) | strncmp('Unknown_',info(1,:),8);
info(:,ind) = [];
sn = 73400000 + round(rand(1)*10^5);
p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'Modality','OT',@ischar);
addParameter(p,'ImageType','ORIGINAL\PRIMARY\AXIAL',@ischar);
addParameter(p,'PatientPosition','HFS',@ischar);
addParameter(p,'PatientID','Unknown',@ischar);
addParameter(p,'StudyDescription','',@ischar);
addParameter(p,'SeriesDescription','',@ischar);
addParameter(p,'StudyID','',@ischar);
addParameter(p,'SeriesNumber',sn,@isnumeric);
addParameter(p,'AcquisitionNumber',1,@isnumeric);
addParameter(p,'ImagePositionPatient',zeros(1,3),@(x)isnumeric(x)&&(numel(x)==3));
addParameter(p,'ImageOrientationPatient',[0;1;0;1;0;0],@(x)isnumeric(x)&&(numel(x)==6));
addParameter(p,'PixelSpacing',ones(1,2),@(x)isnumeric(x)&&(numel(x)==2));
addParameter(p,'SpacingBetweenSlices',1,@isvector);
addParameter(p,'SliceThickness',1,@isvector);
addParameter(p,'RescaleSlope',1,@isvector);
addParameter(p,'RescaleIntercept',0,@isvector);
addParameter(p,'SeriesInstanceUID',dicomuid,@ischar);
addParameter(p,'StudyInstanceUID',dicomuid,@ischar);
parse(p,info{:});
info = p.Results;
info.PatientName = struct('FamilyName',{''},'GivenName',info.PatientID,...
    'MiddleName',{''},'NamePrefix',{''},...
    'NameSuffix',{''});
