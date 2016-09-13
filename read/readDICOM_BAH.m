function [img,label,fov,info] = readDICOM(varargin)
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
    if isfield(info,'PixelSpacing')
        voxsz = info.PixelSpacing(:)';
    else
        voxsz = [1,1];
    end
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
    if isempty(d) || all([nd,nz] == d)
        chk3D = false;
        img = zeros([nd,nz]);
        imgC = cell(1,nz);
        slclocs = zeros(1,nf);
        acqn = zeros(1,nf);
        ithk = zeros(1,nf);
        itpos = zeros(1,nf);
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
            if isfield(info,'ImagePositionPatient')
                slclocs(ifn) = info.ImagePositionPatient(3);
            elseif isfield(info,'SliceLocation')
                slclocs(ifn) = info.SliceLocation;
            elseif isfield(info,'InstanceNumber')
                slclocs(ifn) = info.InstanceNumber;
            else
                slclocs(ifn) = ifn;
            end
            if isfield(info,'AcquisitionNumber') && ~isempty(info.AcquisitionNumber)
                acqn(ifn) = info.AcquisitionNumber;
            end
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

        % Read the file
            timg = double(dicomread(info));
            img(:,:,tframe) = tslp * squeeze(timg) + tint;
            imgC{tframe} = tslp * squeeze(timg) + tint;
        end
        delete(hp);
    else
        disp('Image dimensions do not match')
    end
    if ~isempty(img)
        if chk3D
            
        else
            % Some DICOMs list each slice as a separate acquisition ...
            if (length(unique(acqn(:)))==length(acqn(:))) ...
                    && (length(unique(slclocs(:)))==length(slclocs(:)))
                acqn(:) = acqn(1);
            end
            
            % Separate into acquisitions ...
            iacq = unique(acqn);
            nacq = length(iacq);
            data = struct('iacq',{},'img',{},'slclocs',{},'itpos',{},'ns',{},'thk',{});
            for i = 1:nacq
                data(i).acqn = iacq(i);
                ind = acqn==iacq(i);
                data(i).img = cat(3,imgC{ind});
                data(i).slclocs = slclocs(ind);
                [data(i).slclocs,ord] = sort(data(i).slclocs);
                data(i).img = data(i).img(:,:,ord);
                data(i).itpos = itpos(ind);
                data(i).itpos = data(i).itpos(ord);
                data(i).ns = nnz(ind);
                data(i).thk = ithk(find(ind,1));
            end
            
            % Select acquisition ...
            if nacq>1
                hf = figure;
                nx = ceil(sqrt(nacq));
                ny = ceil(nacq/nx);
                for i = 1:length(data)
                    subplot(ny,nx,i),
                    ri = size(data(i).img,1);
                    hi = imshow(squeeze(data(i).img(round(ri/2),:,:)),[],...
                                'Title',num2str(data(i).iacq));
                    set(get(hi,'Parent'),'DataAspectRatio',[voxsz(2),data(i).thk,1]);
                end
                opts = [{'All (concatenated)'},...
                        cellfun(@strcat,cellfun(@num2str,num2cell(acqID),'UniformOutput',0),...
                                        repmat({':'},1,nacq),...
                                        cellfun(@num2str,num2cell(ns),'UniformOutput',0),...
                                        repmat({' slices'},1,nacq),'UniformOutput',0)];
                answer = listdlg('ListString',opts,'Name','Multiple Acquisitions') - 1;
            end
            
            % Remove redundant slice locations and sort
            A = [acqn(:),slclocs(:),itpos(:)];
            [~,ia,~] = unique(A,'rows');
            if (nf>length(ia)) && ~mod(nf,length(ia))
                % Added 2014-08-14 by BAH
                % MR Solutions .dcm arrayed image data lacks acq # ident
                nv = nf/length(ia);
                for ii = 1:length(ia)
                    ind = ismember(A,A(ia(ii),:),'rows');
                    A(ind,1) = 1:sum(ind);
                end
                [~,ia] = sortrows(A);
            end
            img = img(:,:,ia);
            slclocs = slclocs(ia);
            acqn = acqn(ia);
            ithk = ithk(ia);
            nf = length(ia);
            
            % Reshape to 4D if necessary
            if nf>1
                ns = nf/nv;
                if round(ns)==ns % only if correct # of files loaded
                    img = reshape(img,[nd,ns,nv]);
                else
                    warning('Slices may be off!');
                end
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
        if isfield(info,'SeriesDescription') && ~isempty(info.SeriesDescription)
            label = repmat({info.SeriesDescription},1,nv);
        else
            label = repmat({info.Modality},1,nv);
        end
        if numel(slclocs)>1
            fovz = length(slclocs) * abs(diff(slclocs(1:2)));
        elseif ~isempty(ithk) && (ithk(1)>0)
            thk = ithk(1);
            fovz = nz*thk;
        else
            fovz = nz;
        end
        fov = [double(fliplr(nd)).*voxsz fovz];
        info.SliceVec = slclocs;
    end
end


