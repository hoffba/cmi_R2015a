function [img,label,fov,dcmdata] = readDICOM(varargin)

fnames = varargin{1};
din = [];
if nargin==2
    din = varargin{2};
end
img = []; label = {}; fov = [];
[path,fname,ext] = fileparts(fnames);

% Find and load all DICOM files in selected directory
if strcmpi(fname,'dicomdir')
    % DICOM names are stored in "DICOMDIR" file
    hp = waitbar(0,'Reading DICOMDIR ...');
    tinfo = dicominfo(fnames);
    fnames = fieldnames(tinfo.DirectoryRecordSequence);
    nf = length(fnames);
    ind = false(1,nf);
    for i = 1:nf
        if isfield(tinfo.DirectoryRecordSequence.(fnames{i}),'ReferencedFileID')
            tdir = tinfo.DirectoryRecordSequence.(fnames{i}).ReferencedFileID;
            fnames{i} = strrep(tdir,'\',filesep);
        else
            ind(i) = true;
        end
        waitbar(i/nf,hp,'Finding file locations ...');
    end
    delete(hp);
    fnames(ind) = [];
    nf = nf - nnz(ind);
else
    % All DICOMs are in single directory
    fnames = dir([path,filesep,'*',ext]);
    fnames = {fnames(:).name};
    fnames = fnames(isDICOM(fullfile(path,fnames))==1);
    nf = length(fnames);
end

dcmdata = struct('SeriesInstanceUID',{},... % for sorting multiple images
                 'img',{},...
                 'SlcLoc',{},...
                 'SlcThk',{},...
                 'PixelSpacing',{},...
                 'StudyDescription',{},...
                 'PatientID',{},...
                 'Label',{});
% * Images are sorted into 3D via SlcLoc,
%   and 4D via arrayed values:
%       MRI - TE, TR, ...
%       CT  - kV, ...
             
% Read in all DICOM slices:
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
    
    % Load DICOM file info:
    tinfo = dicominfo(fullfile(path,fnames{ifn}));
    
    % Check if Series already exists in structure:
    j = find(strcmp(tinfo.SeriesInstanceUID,{dcmdata(:).SeriesInstanceUID}),1);
    
    if isempty(j)
        % Initialize new series in structure:
        j = length(dcmdata)+1;
        val = '';
        if isfield(tinfo,'SeriesInstanceUID')
            val = tinfo.SeriesInstanceUID;
        end
        dcmdata(j).SeriesInstanceUID = val;
        val = [1,1];
        if isfield(tinfo,'PixelSpacing')
            val = tinfo.PixelSpacing;
        end
        dcmdata(j).PixelSpacing = val;
        val = [];
        if isfield(tinfo,'SliceThickness')
            val = tinfo.SliceThickness;
        elseif isfield(tinfo,'SpacingBetweenSlices')
            val = tinfo.SpacingBetweenSlices;
        end
        dcmdata(j).SlcThk = val;
        val = '';
        if isfield(tinfo,'StudyDescription')
            val = tinfo.StudyDescription;
        end
        dcmdata(j).StudyDescription = val;
        if isfield(tinfo,'SeriesDescription') && ~isempty(tinfo.SeriesDescription)
            val = tinfo.SeriesDescription;
        elseif isfield(tinfo,'ProtocolName') && ~isempty(tinfo.ProtocolName)
            val = tinfo.ProtocolName;
        else val = tinfo.Modality;
        end
        dcmdata(j).Label = {val};
        dcmdata(j).SlcLoc = nan(nf,1);
        val = '';
        if isfield(tinfo,'PatientID')
            val = tinfo.PatientID;
        end
        dcmdata(j).PatientID = val;
    end
    if isfield(tinfo,'InstanceNumber')
        k = tinfo.InstanceNumber;
    else
        k = length(tinfo.SlcLoc)+1;
    end
    val = k;
    if isfield(tinfo,'SliceLocation')
        val = tinfo.SliceLocation;
    elseif isfield(tinfo,'ImagePositionPatient')
        val = tinfo.ImagePositionPatient;
    end
    dcmdata(j).SlcLoc(k) = val(end);
    if isfield(tinfo,'EchoTime')
        if ~isempty(tinfo.EchoTime)
            dcmdata(j).TE(k) = tinfo.EchoTime;
        end
    end
    if isfield(tinfo,'RepetitionTime')
        if ~isempty(tinfo.RepetitionTime)
            dcmdata(j).TR(k) = tinfo.RepetitionTime;
        end
    end
    
    % Read slice data from file:
    if isfield(tinfo,'RescaleSlope')
        ySlope = tinfo.RescaleSlope;
    else ySlope = 1;
    end
    if isfield(tinfo,'RescaleIntercept')
        yInt = tinfo.RescaleIntercept;
    else yInt = 0;
    end
    timg = ySlope * double(dicomread(tinfo)) + yInt;
    dcmdata(j).img(:,:,k) = timg;
end
delete(hp);

% If multiple series, user selects which to load:
if length(dcmdata)>1
    str = strcat('(',cellfun(@(x)num2str(length(x)),{dcmdata(:).SlcLoc},...
                             'UniformOutput',false),...
                 ' Slices)',[dcmdata(:).Label]);
    answer = listdlg('ListString',str,'SelectionMode','single',...
                     'ListSize',[300,300]);
    if isempty(answer)
        return % User cancelled the load
    else
        dcmdata = dcmdata(answer);
    end
end

ind = isnan(dcmdata.SlcLoc);
dcmdata.img(:,:,ind) = [];
dcmdata.SlcLoc(ind) = [];

% Sort slices by location:
uind = [];
n4d = 1;
if all(isfield(dcmdata,{'TR','TE'}))
    [C,~,uind] = unique([dcmdata.TR',dcmdata.TE'],'rows');
    n4d = size(C,1);
end
uind = [uind,dcmdata.SlcLoc];
[~,ix,~] = unique([uind,dcmdata.SlcLoc],'rows');
ns = length(ix);
[d(1),d(2),d(3)] = size(dcmdata.img);
if ns~=d(3)
    % MR Solutions DICOMs have no way to separate out b-values of DWI
    if mod(d(3),ns)==0
        n4d =  d(3)/ns;
        uind = repmat(1:n4d,ns,1);
        uind = uind(:);
        ix = repmat(ix,n4d,1);
    else
        error('Matrix cannot be separated into 4D.')
    end
end
if (n4d>1) && (d(3)==(n4d*nnz(uind==1)))
    ns = d(3)/n4d;
    dcmdata.img = reshape(dcmdata.img(:,:,ix),[d(1:2),ns,n4d]);
    dcmdata.SlcLoc = dcmdata.SlcLoc(ix(1:ns));
    % Determine labeling:
    dcmdata.Label = cell(1,n4d);
    ind = ix(1:ns:d(3));
    if length(unique(dcmdata.TR))>1
        dcmdata.Label = strcat(dcmdata.Label,'TR=',...
            cellfun(@num2str,num2cell(dcmdata.TR(ind)),...
                    'UniformOutput',false),'; ');
    end
    if length(unique(dcmdata.TE))>1
        dcmdata.Label = strcat(dcmdata.Label,'TE=',...
            cellfun(@num2str,num2cell(dcmdata.TE(ind)),...
                    'UniformOutput',false),'; ');
    end
    ind = cellfun(@isempty,dcmdata.Label);
    if any(ind)
        dcmdata.Label(ind) = cellfun(@num2str,num2cell(find(ind)));
        dcmdata.Label = cellfun(@num2str,num2cell(1:n4d),...
                                'UniformOutput',false);
    end
end

% if all(n==n(1)) % (all slice locations replicated the same number of times)
%         dcmdata.img = reshape(dcmdata.img,[d(1:2),d(3)/n(1),n(1)]);
% end
% if isfield(tinfo,'SequenceName') && strncmp(tinfo.SequenceName,'MGE',3)
%     nz = nnz(dcmdata.SlcLoc==dcmdata.SlcLoc(1));
%     ns = d(3)/nz;
%     dcmdata.img = reshape(dcmdata.img(:,:,ix),[d(1:2),ns,nz]);
% else
%     % MR Solutions DICOMs aren't separated with AcquisitionNumber:
%     n = zeros(1,length(ix));
%     for i = 1:length(ix)
%         n(i) = nnz(ic==ix(i));
%     end
%     if all(n==n(1))
%         dcmdata.img = reshape(dcmdata.img,[d(1:2),d(3)/n(1),n(1)]);
%         dcmdata.Label = repmat(dcmdata.Label,1,n(1));
%     end
%     dcmdata.img = dcmdata.img(:,:,ix,:);
%     dcmdata.SlcLoc = dcmdata.SlcLoc(ix);
% end

% Check for dimension match
[d(1),d(2),d(3),d(4)] = size(dcmdata.img);
if ~isempty(din) && ~all(d(1:3)==din)
    error(['Dimensions must match!',sprintf(' %f',d)])
end

% Prepare data to return:
img = dcmdata.img;
dcmdata = rmfield(dcmdata,'img');
[d(1),d(2),d(3),d(4)] = size(img);
fov = [fliplr(dcmdata.PixelSpacing(:)'),1] .* d(1:3);
if length(dcmdata.SlcLoc)>1
    fov(3) = abs(diff(dcmdata.SlcLoc([1,end]))) * d(3)/(d(3)-1);
elseif ~isempty(dcmdata.SlcThk)
    fov(3) = d(3) * dcmdata.SlcThk;
end
label = dcmdata.Label;

% % Read first DICOM for header info
% info = dicominfo(fullfile(path,fnames{1}));
% % Needed header fields:
% nfields = {'Rows','Columns','PixelSpacing','Modality'};
% if all(isfield(info,nfields))
%     nd = [info.Rows,info.Columns];
%     if isfield(info,'PixelSpacing')
%         voxsz = info.PixelSpacing(:)';
% %     elseif isfield(info,
%     else
%         voxsz = [1,1];
%     end
%     if all(isfield(info,{'NumberOfTemporalPositions',...
%                          'TemporalPositionIdentifier'}))
%         nv = info.NumberOfTemporalPositions;
%     else
%         nv = 1;
%     end
%     if isfield(info,'NumberOfFrames')
%         nz = info.NumberOfFrames;
%     else
%         nz = nf;
%     end
%     if isempty(d) || all([nd,nz] == d)
%         chk3D = false;
%         img = zeros([nd,nz]);
%         slclocs = zeros(1,nf);
%         acqn = zeros(1,nf);
%         ithk = zeros(1,nf);
%         itpos = zeros(1,nf);
%         hp = waitbar(0,'','WindowStyle','modal',...
%             'CreateCancelBtn','setappdata(gcbf,''canceling'',1)'); % Option to cancel load midway
%         setappdata(hp,'canceling',0)
%         for ifn = 1:nf
%         % Check for Cancel button press
%             if getappdata(hp,'canceling')
%                 img = []; label = {}; fov = [];
%                 break
%             end
%             waitbar((ifn-1)/nf,hp,fnames{ifn});
%             info = dicominfo(fullfile(path,fnames{ifn}));
%         % Check for rescale slope and intercept:
%             tslp = 1; tint = 0; % Defaults if fields not found
%             if isfield(info,'RescaleSlope')
%                 tslp = info.RescaleSlope;
%             end
%             if isfield(info,'RescaleIntercept')
%                 tint = info.RescaleIntercept;
%             end
%             if isfield(info,'ImagePositionPatient')
%                 slclocs(ifn) = info.ImagePositionPatient(3);
%             elseif isfield(info,'SliceLocation')
%                 slclocs(ifn) = info.SliceLocation;
%             elseif isfield(info,'InstanceNumber')
%                 slclocs(ifn) = info.InstanceNumber;
%             else
%                 slclocs(ifn) = ifn;
%             end
%             if isfield(info,'AcquisitionNumber') && ~isempty(info.AcquisitionNumber)
%                 acqn(ifn) = info.AcquisitionNumber;
%             end
%             if isfield(info,'SliceThickness') && ~isempty(info.SliceThickness)
%                 ithk(ifn) = info.SliceThickness;
%             elseif isfield(info,'SpacingBetweenSlices')
%                 ithk(ifn) = info.SpacingBetweenSlices;
%             end
%             if isfield(info,'TemporalPositionIdentifier')
%                 itpos(ifn) = info.TemporalPositionIdentifier;
%             end
%             if isfield(info,'NumberOfFrames') % 3D image
%                 tframe = 1:info.NumberOfFrames;
%                 chk3D = true;
%             else
%                 tframe = ifn;
%             end
% 
%         % Read the file
%             if (info.Rows==nd(1)) && (info.Columns==nd(2))
%                 timg = double(dicomread(info));
%                 img(:,:,tframe) = tslp * squeeze(timg) + tint;
%             end
%         end
%         delete(hp);
%     else
%         disp('Image dimensions do not match')
%     end
%     if ~isempty(img)
%         if chk3D
%             
%         else
%             % Some DICOMs list each slice as a separate acquisition ...
%             if (length(unique(acqn(:)))==length(acqn(:))) ...
%                     && (length(unique(slclocs(:)))==length(slclocs(:)))
%                 acqn(:) = acqn(1);
%             end
%             
%             % Remove redundant slice locations and sort
%             A = [acqn(:),slclocs(:),itpos(:)];
%             [~,ia,~] = unique(A,'rows');
%             if (nf>length(ia)) && ~mod(nf,length(ia))
%                 % Added 2014-08-14 by BAH
%                 % MR Solutions .dcm arrayed image data lacks acq # ident
%                 nv = nf/length(ia);
%                 uacqn = unique(acqn);
%                 inacq = zeros(1,length(uacqn));
%                 for ii = 1:length(uacqn)
%                     inacq(ii) = sum(acqn==uacqn(ii));
%                 end
%                 if all(inacq==inacq(1))
%                     for ii = 1:length(ia)
%                         ind = ismember(A,A(ia(ii),:),'rows');
%                         A(ind,1) = 1:sum(ind);
%                     end
%                     [~,ia] = sortrows(A);
%                 else % Found some redundant images in Van Poznak's data
%                     nv = 1;
%                 end
%             end
%             img = img(:,:,ia);
%             slclocs = slclocs(ia);
%             acqn = acqn(ia);
%             ithk = ithk(ia);
%             nf = length(ia);
%             
%             % Reshape to 4D if necessary
%             if nf>1
%                 ns = nf/nv;
%                 if round(ns)==ns % only if correct # of files loaded
%                     img = reshape(img,[nd,ns,nv]);
%                 else
%                     warning('Slices may be off!');
%                 end
%             end
%         end
%         
%         % Choose which acquisition to load:
%         acqID = unique(acqn); nacq = length(acqID);
%         if nacq>1
%             hf = figure;
%             disp(num2str(slclocs(:)));
%             imin = min(img(:)); imax = max(img(:));
%             imshow((squeeze(img(round(size(img,1)/2),:,:,1))-imin)/(imax-imin));
%             ns = zeros(1,nacq);
%             for i = 1:nacq
%                 ns(i) = sum(acqID(i)==acqn);
%             end
%             opts = [{'All (concatenated)'},...
%                     cellfun(@strcat,cellfun(@num2str,num2cell(acqID),'UniformOutput',0),...
%                                     repmat({':'},1,nacq),...
%                                     cellfun(@num2str,num2cell(ns),'UniformOutput',0),...
%                                     repmat({' slices'},1,nacq),'UniformOutput',0)];
%             answer = listdlg('ListString',opts,'Name','Multiple Acquisitions') - 1;
%             close(hf);
%         else
%             answer = 1;
%         end
%         % Remove unwanted slices:
%         if ~isempty(answer) && all(answer>0)
%             ind = ~ismember(acqn,acqID(answer));
%             img(:,:,ind,:) = [];
%             slclocs(ind) = [];
%             ithk(ind) = [];
%         end
%         if isfield(info,'SeriesDescription') && ~isempty(info.SeriesDescription)
%             label = repmat({info.SeriesDescription},1,nv);
%         else
%             label = repmat({info.Modality},1,nv);
%         end
%         if numel(slclocs)>1
%             fovz = length(slclocs) * abs(diff(slclocs(1:2)));
%         elseif ~isempty(ithk) && (ithk(1)>0)
%             thk = ithk(1);
%             fovz = nz*thk;
%         else
%             fovz = nz;
%         end
%         fov = [double(fliplr(nd)).*voxsz fovz];
%         info.SliceVec = slclocs;
%     end
% end


