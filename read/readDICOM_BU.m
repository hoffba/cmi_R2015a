function [img,label,fov] = readDICOM(varargin)
fnames = varargin{1};
[path,~,ext] = fileparts(fnames);
% Find and load all DICOM files in selected directory
fnames = dir([path filesep '*' ext]);
if isempty(ext)
    fnames(1:2) = [];
end
nf = length(fnames);
% Read first DICOM for header info
info = dicominfo([path filesep fnames(1).name]);
% Make sure DICOM has all needed header attributes
if all(isfield(info,{'Rows','Columns','SliceLocation','Modality','PixelSpacing'}))
    % Lung CT
    nd = [info.Rows,info.Columns];
    img = zeros([nd nf]);
    slclocs = zeros(1,nf);
    tfind = zeros(1,nf);
    hp = waitbar(0,'','WindowStyle','modal',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)'); % Option to cancel load midway
    setappdata(hp,'canceling',0)
    tic
    for ifn = 1:nf
        % Check for Cancel button press
        if getappdata(hp,'canceling')
            img = []; label = {}; fov = [];
            break
        end
        waitbar((ifn-1)/nf,hp,fnames(ifn).name);
        fpath = [path filesep fnames(ifn).name];
        info = dicominfo(fpath);
        tslp = 1; tint = 0; % Defaults if fields not found
        if isfield(info,'RescaleSlope')
            tslp = info.RescaleSlope;
        end
        if isfield(info,'RescaleIntercept')
            tint = info.RescaleIntercept;
        end
        slclocs(ifn) = info.SliceLocation;
        tfind(ifn) = info.AcquisitionNumber;
        img(:,:,ifn) = tslp * double(dicomread(info)) + tint;
    end
    delete(hp);
    if ~isempty(img)
        [slclocs,order] = sort(slclocs,'descend');
        img = img(:,:,fliplr(order));
        toc
        label = {info.Modality};
        fov = [double(nd).*info.PixelSpacing' 2*slclocs(1)-slclocs(end)-slclocs(2)];
    end
elseif all(isfield(info,{'Rows','Columns','PixelSpacing','Modality'})) && isempty(info.Modality)
    % ImBio segmentation DICOM (mask)
    % Assumes three regions: Trachea, Left Lung, Right Lung
    % Removes trachea region from the image
    nd = [info.Rows,info.Columns];
    img = zeros([nd nf]);
    hp = waitbar(0,'','WindowStyle','modal',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)'); % Option to cancel load midway
    setappdata(hp,'canceling',0)
    tic
    for ifn = 1:nf
        % Check for Cancel button press
        if getappdata(hp,'canceling')
            img = []; label = {}; fov = [];
            break
        end
        waitbar((ifn-1)/nf,hp,fnames(ifn).name);
        fpath = [path filesep fnames(ifn).name];
        info = dicominfo(fpath);
        img(:,:,ifn) = double(dicomread(info)) * info.RescaleSlope ...
                                               + info.RescaleIntercept;
    end
    delete(hp);
    if ~isempty(img)
        toc
        label = {'ImBio-Seg'};
        thk = 0.75; % Assume a slice thickness for the data since none is given
        fov = [double(nd).*info.PixelSpacing' thk*nf];
    end
elseif all(isfield(info,{'Rows','Columns','PixelSpacing','SliceThickness','ImageIndex','Modality',...
        'RescaleIntercept','RescaleSlope'}))
    % Siemens Inveon Data
    nd = [info.Rows,info.Columns];
    img = zeros([nd nf]);
    hp = waitbar(0,'','WindowStyle','modal',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)'); % Option to cancel load midway
    setappdata(hp,'canceling',0)
    tic
    for ifn = 1:nf
        % Check for Cancel button press
        if getappdata(hp,'canceling')
            img = []; label = {}; fov = [];
            break
        end
        waitbar((ifn-1)/nf,hp,fnames(ifn).name);
        fpath = [path filesep fnames(ifn).name];
        info = dicominfo(fpath);
        img(:,:,info.ImageIndex) = double(dicomread(info))*info.RescaleSlope + info.RescaleIntercept;
    end
    delete(hp);
    if ~isempty(img)
        toc
        label = {info.Modality};
        fov = [double(nd).*info.PixelSpacing' info.SliceThickness*nf];
    end
else
    error('DICOM header does not have required attributes!')
end
