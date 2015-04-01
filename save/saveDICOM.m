function status = saveDICOM(fname,img,label,max_ext,info,verbchk)
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
    
    % Determine DICOM header tags:
    if (nargin<5)
        info = [];
    end
    sn = 73400000 + round(rand(1)*10^5);
    if ~isfield(info,'PatientID')
        info.PatientID = fname;
    end
    if ~isfield(info,'Modality')
        info.Modality = 'CT';
    end
    if ~isfield(info,'StudyDescription')
        info.StudyDescription = 'Matlab Analysis';
    end
    if ~isfield(info,'StudyID')
        info.StudyID = num2str(sn);
    end
    if ~isfield(info,'SeriesNumber')
        info.SeriesNumber = sn;
    end
    if ~isfield(info,'AcquisitionNumber')
        info.AcquisitionNumber = 1;
    end
    if ~isfield(info,'ImagePositionPatient')
        info.ImagePositionPatient = zeros(1,3);
    end
    if ~isfield(info,'ImageOrientationPatient')
        info.ImageOrientationPatient = [0,1,0,1,0,0];
    end
    
    rtags = {'PatientID','Modality','StudyID','StudyDescription'};
    if (nargin<5) || ~all(isfield(info,rtags))
        prompt = {'Modality (MR,CT,PT=PET,ST=SPECT,SEG=Segmentation,OT=other)',...
                  'PatientID:',...
                  'StudyID',...
                  'StudyDescription:'};
        answer = {info.Modality,info.PatientID,info.StudyID,info.StudyDescription};
        answer = inputdlg(prompt,'DICOM Tags',1,answer);
    else
        answer = true;
    end
    
    if ~isempty(answer)
        if iscellstr(answer)
            info.Modality = answer{1};
            info.PatientID = answer{2};
            info.StudyID = answer{3};
            info.StudyDescription = answer{4};
        end
    
        % Auto-scale to int16:
        mxscale = 2^15-1;
        sthk = max_ext(3)/ns;
        spos = (sthk - max_ext(3))/2;
        warning('off','images:dicomwrite:inconsistentIODAndCreateModeOptions');
        
        % Initialize DICOM info:
        info.PatientName = struct('FamilyName',{''},'GivenName',info.PatientID,...
                                'MiddleName',{''},'NamePrefix',{''},...
                                'NameSuffix',{''});
        info.SeriesDescription = '';
        info.PixelSpacing = max_ext(1:2)./d;
        info.SpacingBetweenSlices = sthk;
        info.SliceThickness = sthk;
        info.ImageOrientationPatient = [0,1,0,1,0,0];
        info.RescaleSlope = 1;
        info.RescaleIntercept = 0;

%         if (nargin<5) && verbchk
%             answer = questdlg('Save as 2D or 3D?','DICOM save','2D','3D','Cancel','3D');
%         else
            answer = '2D';
%         end
        for ivec = 1:nv
            timg = permute(img(:,:,:,ivec),[1,2,4,3]);
            tmax = max(timg(:));
            tmin = min(timg(:));
            info.RescaleIntercept = tmin;
            info.RescaleSlope = (tmax - tmin)/mxscale;
            if info.RescaleSlope == 0
                info.RescaleSlope = 1;
            end
            timg = int16((timg-info.RescaleIntercept)/info.RescaleSlope);
            info.SeriesDescription = label{ivec};
            info.ImagePositionPatient = [0,0,spos];
            info.SeriesInstanceUID = dicomuid;
            info.StudyInstanceUID = dicomuid;

            if nv>1
                fnout = [fname,'_',num2str(ivec)];
            else
                fnout = fname;
            end

            if strcmp(answer,'2D')
                if ~exist(fullfile(pathstr,fnout),'dir')
                    mkdir(fullfile(pathstr,fnout));
                end
                if verbchk, hw = waitbar(0,'Saving DICOM slices ...'); end
                for islc = 1:ns
                    info.ImagePositionPatient(3) = info.ImagePositionPatient(3) + sthk;
                    info.SliceLocation = info.ImagePositionPatient(3);
                    info.InstanceNumber = islc;
                    stat = dicomwrite(timg(:,:,1,islc),...
                                fullfile(pathstr,fnout,sprintf([fnout,'_s%04u.dcm'],islc)),...
                                info);
                    status = status && isempty(stat.BadAttribute) ...
                                    && isempty(stat.MissingCondition) ...
                                    && isempty(stat.MissingData) ...
                                    && isempty(stat.SuspectAttribute);
                    if verbchk, waitbar(islc/ns,hw,['Slc ',num2str(info.ImagePositionPatient(3))]); end
                end
                if verbchk, close(hw); end
            elseif strcmp(answer,'3D')
                stat = dicomwrite(timg,...
                            fullfile(pathstr,[fnout,'.dcm']),info,...
                            'MultiframeSingleFile',true);
                status = status && isempty(stat.BadAttribute) ...
                                && isempty(stat.MissingCondition) ...
                                && isempty(stat.MissingData) ...
                                && isempty(stat.SuspectAttribute);
            end
        end
    end
end


