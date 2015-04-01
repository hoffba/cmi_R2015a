function [img,label,fov] = readDICOM(varargin)
fnames = varargin{1};
d = varargin{2};
img = []; label = {}; fov = [];
[path,~,ext] = fileparts(fnames);
% Find and load all DICOM files in selected directory
fnames = dir([path,filesep,'*',ext]);
fnames = {fnames(:).name};
fnames = fnames(isDICOM(fullfile(path,fnames))==1);
% if isempty(ext)
%     fnames2 = dir([path,filesep,'*.*']);
%     fnames2 = {fnames2(:).name};
%     fnames = setdiff(fnames,fnames2);
% end
nf = length(fnames);
% Read first DICOM for header info
info = dicominfo(fullfile(path,fnames{1}));
% Needed header fields:
nfields = {'Rows','Columns','PixelSpacing','Modality'};
if all(isfield(info,nfields))
    nd = [info.Rows,info.Columns];
    if all(isfield(info,{'NumberOfTemporalPositions',...
                         'TemporalPositionIdentifier'}))
        nv = info.NumberOfTemporalPositions;
    else
        nv = 1;
    end
    if isfield(info,'NumberOfFrames')
        nz = info.NumberOfFrames;
    else
        nz = nf;
    end
    if isempty(d) || all([nd,nf] == d)
        chk3D = false;
        img = zeros([nd,nz]);
        slclocs = zeros(1,nf);
        acqn = zeros(1,nf);
        ithk = zeros(1,nf);
%         ind = false(1,nf);
        itpos = zeros(1,nf);
    %tinfo(nf) = info;
        hp = waitbar(0,'','WindowStyle','modal',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)'); % Option to cancel load midway
        setappdata(hp,'canceling',0)
        for ifn = 1:nf
        % Check for Cancel button press
            if getappdata(hp,'canceling')
                img = []; label = {}; fov = [];
                break
            end
            waitbar((ifn-1)/nf,hp,fnames{ifn});
            info = dicominfo(fullfile(path,fnames{ifn}));
        % Check for rescale slope and intercept:
            tslp = 1; tint = 0; % Defaults if fields not found
            if isfield(info,'RescaleSlope')
                tslp = info.RescaleSlope;
            end
            if isfield(info,'RescaleIntercept')
                tint = info.RescaleIntercept;
            end
            if isfield(info,'SliceLocation')
                slclocs(ifn) = info.SliceLocation;
            elseif isfield(info,'InstanceNumber')
                slclocs(ifn) = info.InstanceNumber;
            elseif isfield(info,'ImagePositionPatient')
                slclocs(ifn) = info.ImagePositionPatient(3);
            else
                slclocs(ifn) = ifn;
            end
            if isfield(info,'AcquisitionNumber') && ~isempty(info.AcquisitionNumber)
                acqn(ifn) = info.AcquisitionNumber;
            end
%             if isfield(info,'OverlayData_0')
%                 ind(ifn) = true;
%             end
            if isfield(info,'SliceThickness')
                ithk(ifn) = info.SliceThickness;
            elseif isfield(info,'SpacingBetweenSlices')
                ithk(ifn) = info.SpacingBetweenSlices;
            end
            if isfield(info,'TemporalPositionIdentifier')
                itpos(ifn) = info.TemporalPositionIdentifier;
            end
            if isfield(info,'NumberOfFrames') % 3D image
                tframe = 1:info.NumberOfFrames;
                chk3D = true;
            else
                tframe = ifn;
            end
    %tinfo(ifn) = info;

        % Read the file
            timg = double(dicomread(info));
%             if isfield(info,'PhotometricInterpretation') && ...
%                     strcmp(info.PhotometricInterpretation,'RGB')
%                 timg = rgb2ind(timg(:,:,1:3),2^16);
%             end
            img(:,:,tframe) = tslp * squeeze(timg) + tint;
        end
        delete(hp);
    else
        disp('Image dimensions do not match')
    end
    if ~isempty(img)
        if chk3D
            
        else
%             % Remove overlay images
%             nover = sum(ind); disp(['# overlays = ' num2str(nover)])
%             if any(ind)
%                 img(:,:,ind,:) = [];
%                 slclocs(ind) = [];
%                 acqn(ind) = [];
%                 ithk(ind) = [];
%             end
            % Remove redundant slice locations and sort
            [~,ia,~] = unique([acqn(:),itpos(:),slclocs(:)],'rows');
%             img = img(:,:,ia,:);
%             slclocs = slclocs(ia);
%             acqn = acqn(ia);
%             ithk = ithk(ia);
%             nf = length(ia);
%             % Sort slices into position:
%             [~,order] = sort(slclocs);
%             % Sort into acquisition order
%             [acqn,torder] = sort(acqn(order));
%             order = order(torder);
            % Reshape to 4D if necessary
            if nf>1
                ns = nf/nv;
                if round(ns)==(ns) % only if correct # of files loaded
                    [~,torder] = sort(itpos(ia));
                    ia = ia(torder);
                    img = reshape(img(:,:,ia),[nd,ns,nv]);
                else
                    ns = nf; % Don't separate into 4D
                    img = img(:,:,ia);
                end
                % Adjust data into order
                ithk = ithk(ia(1:ns));
                slclocs = slclocs(ia(1:ns));
            end
        end
        % Choose which acquisition to load:
        acqID = unique(acqn); nacq = length(acqID);
        if nacq>1
            hf = figure;
            disp(num2str(slclocs(:)));
            imin = min(img(:)); imax = max(img(:));
            imshow((squeeze(img(round(size(img,1)/2),:,:,1))-imin)/(imax-imin));
            ns = zeros(1,nacq);
            for i = 1:nacq
                ns(i) = sum(acqID(i)==acqn);
            end
            opts = [{'All (concatenated)'},...
                    cellfun(@strcat,cellfun(@num2str,num2cell(acqID),'UniformOutput',0),...
                                    repmat({':'},1,nacq),...
                                    cellfun(@num2str,num2cell(ns),'UniformOutput',0),...
                                    repmat({' slices'},1,nacq),'UniformOutput',0)];
            answer = listdlg('ListString',opts,'Name','Multiple Acquisitions') - 1;
            close(hf);
        else
            answer = 1;
        end
        % Remove unwanted slices:
        if ~isempty(answer) && all(answer>0)
            ind = ~ismember(acqn,acqID(answer));
            img(:,:,ind,:) = [];
            slclocs(ind) = [];
            ithk(ind) = [];
        end
        label = repmat({info.Modality},1,nv);
        if numel(slclocs)>1
            thk = slclocs(2)-slclocs(1);
            fovz = slclocs(end)-slclocs(1)+thk;
        elseif ~isempty(ithk) && (ithk(1)>0)
            thk = ithk(1);
            fovz = nz*thk;
        else
            fovz = nz;
        end
        fov = [double(fliplr(nd)).*info.PixelSpacing' fovz];
    end
end


