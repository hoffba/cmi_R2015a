function status = saveDICOMnowaitbar(fname,img,label,max_ext,info)
%% Save as 3D DICOM
status = false;
if (nargin>=4)
    
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
    rtags = {'PatientID','Modality','StudyDescription'};
    if (nargin==5) && all(isfield(info,rtags))
        status = true;
    end
    
    if status
        % Auto-scale to int16:
        mxscale = 2^15-1;
        spos = max_ext./(2*[d,ns]);
        sthk = max_ext(3)/ns;
        warning('off','images:dicomwrite:inconsistentIODAndCreateModeOptions');
        info = struct('PatientID',info.PatientID,...
                      'PatientName',struct('FamilyName',info.PatientID),...
                      'StudyDescription',info.StudyDescription,...
                      'SeriesDescription','',...
                      'PixelSpacing',max_ext(1:2)./d,...
                      'SpacingBetweenSlices',max_ext(3)/ns,...
                      'SliceThickness',sthk,...
                      'ImageOrientationPatient',[0,1,0,1,0,0],...
                      'Modality',info.Modality,...
                      'RescaleSlope',1,...
                      'RescaleIntercept',0);

        answer = '2D';
        for ivec = 1:nv
            timg = permute(img(:,:,:,ivec),[1,2,4,3]);
            tmax = max(timg(:));
            tmin = min(timg(:));
            info.RescaleIntercept = tmin;
            info.RescaleSlope = (tmax - tmin)/mxscale;
            timg = int16((timg-info.RescaleIntercept)/info.RescaleSlope);
            info.SeriesDescription = label{ivec};
            info.ImagePositionPatient = -spos;
            info.SeriesInstanceUID = dicomuid;
            info.StudyInstanceUID = dicomuid;

            fnout = [fname,'_',num2str(ivec)];

            if strcmp(answer,'2D')
                if ~exist(fullfile(pathstr,fnout),'dir')
                    mkdir(fullfile(pathstr,fnout));
                end
                for islc = 1:ns
                    info.ImagePositionPatient(3) = info.ImagePositionPatient(3) + sthk;
                    info.SliceLocation = (islc-1)*info.SpacingBetweenSlices;
                    info.InstanceNumber = islc;
                    stat = dicomwrite(timg(:,:,1,islc),...
                                fullfile(pathstr,fnout,sprintf([fnout,'_s%04u.dcm'],islc)),...
                                info);
                    status = status && isempty(stat.BadAttribute) ...
                                    && isempty(stat.MissingCondition) ...
                                    && isempty(stat.MissingData) ...
                                    && isempty(stat.SuspectAttribute);
                end
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


