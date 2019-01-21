function [img,label,fov,dcmdata] = readDICOM(varargin)
% Assumes rectilinear geometry

fnames = varargin{1};
din = [];
if nargin==2
    din = varargin{2};
end
img = []; label = {}; fov = [];
[fpath,fname,ext] = fileparts(fnames);
if isempty(fpath)
    fpath = pwd;
end
if ~any(strcmp(ext,{'.dcm','.1'}))
    ext = [];
end
if strncmp(version('-release'),'2016',4)
    v = {'UseDictionaryVR',true};
else
    v = {};
end

% Find and load all DICOM files in selected directory
if strcmpi(fname,'dicomdir')
    % DICOM names are stored in "DICOMDIR" file
    hp = waitbar(0,'Reading DICOMDIR ...');
    tinfo = dicominfo(fnames,v{:});
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
    fnames = dir([fpath,filesep,'*',ext]);
    fnames = {fnames(:).name};
    fnames(~cellfun(@(x)(exist(fullfile(fpath,x),'file')==2)&&isdicom(fullfile(fpath,x)),fnames)) = [];
%     fnames = fnames(isDICOM(fullfile(path,fnames))==1);
    nf = length(fnames);
end

dcmdata = struct('SeriesInstanceUID',{},... % for sorting multiple images
                 'img',{},...
                 'd',{},...
                 'AcquisitionNumber',{},...
                 'SeriesNumber',{},...
                 'TemporalPositionIdentifier',{},...
                 'DiffusionNumber',{},...
                 'PixelSpacing',{},...
                 'SliceOrient',{},...
                 'SlicePos',{},...
                 'SlcThk',{},...
                 'StudyDescription',{},...
                 'PatientID',{},...
                 'Label',{});
% * Images are sorted into 3D via SlcLoc,
%   and 4D via arrayed values:
%       MRI - TE, TR, ...
%       CT  - kV, ...
             
% Read in all DICOM slices:
hp = waitbar(0,'','Name','Loading DICOM image ...','WindowStyle','modal',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)'); % Option to cancel load midway
setappdata(hp,'canceling',0)
set(hp.Children(2).Title,'Interpreter','none');
try
for ifn = 1:nf
    
    % Check for Cancel button press
    if getappdata(hp,'canceling')
        img = []; label = {}; fov = [];
        dcmdata = [];
        break
    end
    waitbar((ifn-1)/nf,hp,fnames{ifn});
    
    % Load DICOM file info:
    tinfo = dicominfo(fullfile(fpath,fnames{ifn}),v{:});
    if (ifn==1) % Check manufacturer to read private tags:
        ddir = fileparts(which('cmi'));
        if isfield(tinfo,'Manufacturer')
            str = tinfo.Manufacturer;
        else
            str = '';
        end
        if strncmpi(str,'philips',7)
            dicomdict('set',fullfile(ddir,'read','dicom-dict-Philips.txt'));
        elseif strncmpi(str,'ge',2)
            dicomdict('set',fullfile(ddir,'read','dicom-dict-GE.txt'));
        elseif strncmpi(str,'siemens',7)
            dicomdict('set',fullfile(ddir,'read','dicom-dict-Siemens.txt'));
        else
            dicomdict('factory');
        end
        tinfo = dicominfo(fullfile(fpath,fnames{ifn}),v{:});
    end
    
    if isfield(tinfo,'Modality')
        % Check if Series already exists in structure:
    %     disp(tinfo.SeriesInstanceUID);
    %     disp(num2str(tinfo.AcquisitionNumber));
    %     disp(num2str(tinfo.SeriesNumber));
    %     disp(num2str(tinfo.InstanceNumber));
    %     waitforbuttonpress
%         acqN = 0;
%         if isfield(tinfo,'AcquisitionNumber') && ~isempty(tinfo.AcquisitionNumber)
%             acqN = tinfo.AcquisitionNumber;
%         end
        serN = 0;
        if isfield(tinfo,'SeriesNumber') && ~isempty(tinfo.SeriesNumber)
            serN = tinfo.SeriesNumber;
        end
        tempN = 0;
        if isfield(tinfo,'TemporalPositionIdentifier') && ~isempty(tinfo.TemporalPositionIdentifier)
            tempN = tinfo.TemporalPositionIdentifier;
        end
        diffN = 0;
        if isfield(tinfo,'DiffusionBValue')
            nv = length(tinfo.DiffusionBValue);
            if nv==1
                diffN = tinfo.DiffusionBValue;
            elseif isa(tinfo.DiffusionBValue,'uint8')
                diffN = str2double(char(tinfo.DiffusionBValue'));
            end
        elseif isfield(tinfo,'Private_DiffusionBValue')
            nv = length(tinfo.Private_DiffusionBValue);
            if nv==1
                diffN = tinfo.Private_DiffusionBValue;
            elseif isa(tinfo.Private_DiffusionBValue,'uint8')
                diffN = str2double(char(tinfo.Private_DiffusionBValue'));
            end
        elseif isfield(tinfo,'B_value') && (length(tinfo.B_value)==1)
            diffN = tinfo.B_value;
        elseif isfield(tinfo,'SlopInt_6_9') && ~isempty(tinfo.SlopInt_6_9)
            if isa(tinfo.SlopInt_6_9,'uint8')
                % stored as char
                diffN = char(tinfo.SlopInt_6_9');
                diffN = str2double(strtok(diffN,'\'));
                if diffN>1e6
                    diffN = diffN - 1e9;
                end
            else % assumed double
                diffN = tinfo.SlopInt_6_9(1);
            end
        end
        diffD = zeros(3,1);
        if isfield(tinfo,'DiffusionGradientOrientation')
            nv = length(tinfo.DiffusionGradientOrientation);
            if nv==3
                diffD = tinfo.DiffusionGradientOrientation;
            elseif isa(tinfo.DiffusionGradientOrientation,'uint8') && (nv>=24)
                diffD = [ typecast(tinfo.DiffusionGradientOrientation(1:8),'double') ;...
                          typecast(tinfo.DiffusionGradientOrientation(9:16),'double') ;...
                          typecast(tinfo.DiffusionGradientOrientation(17:24),'double') ];
            end
        elseif isfield(tinfo,'DiffusionGradientDirection') && (length(tinfo.DiffusionGradientDirection)==3)
            diffD = tinfo.DiffusionGradientDirection;
        elseif all(isfield(tinfo,{'DiffusionDirectionX','DiffusionDirectionY','DiffusionDirectionZ'}))
            if ~isempty(tinfo.DiffusionDirectionX) && (numel(tinfo.DiffusionDirectionX)==1)
                diffD(1) = tinfo.DiffusionDirectionX;
            end
            if ~isempty(tinfo.DiffusionDirectionY) && (numel(tinfo.DiffusionDirectionY)==1)
                diffD(2) = tinfo.DiffusionDirectionY;
            end
            if ~isempty(tinfo.DiffusionDirectionZ) && (numel(tinfo.DiffusionDirectionZ)==1)
                diffD(3) = tinfo.DiffusionDirectionZ;
            end
        end
        TE = 0;
        if isfield(tinfo,'EchoTime') && ~isempty(tinfo.EchoTime)
            TE = tinfo.EchoTime;
        end
        TR = 0;
        if isfield(tinfo,'RepetitionTime') && ~isempty(tinfo.RepetitionTime)
            TR = tinfo.RepetitionTime;
        end
        j = [];
        if ~isempty(dcmdata)
            j = find( strcmp(tinfo.SeriesInstanceUID,{dcmdata(:).SeriesInstanceUID}) ...
                    & (serN==[dcmdata(:).SeriesNumber]) ...
                    & (tempN == [dcmdata(:).TemporalPositionIdentifier]) ...
                    & (diffN == [dcmdata(:).DiffusionNumber]) ...
                    & ismember([dcmdata(:).DiffusionDir]',diffD','rows')' ...
                    & (TE==[dcmdata(:).TE]) ...
                    & (TR==[dcmdata(:).TR]) ,1);
        end
        if isempty(j)
            % Initialize new series in structure:
            j = length(dcmdata)+1;
            val = '';
            if isfield(tinfo,'SeriesInstanceUID')
                val = tinfo.SeriesInstanceUID;
            end
            dcmdata(j).SeriesInstanceUID = val;
            dcmdata(j).SeriesNumber = serN;
            dcmdata(j).TemporalPositionIdentifier = tempN;
            dcmdata(j).DiffusionNumber = diffN;
            dcmdata(j).DiffusionDir = diffD;
            dcmdata(j).TE = TE;
            dcmdata(j).TR = TR;
            val = [1,1];
            if isfield(tinfo,'PixelSpacing')
                val = tinfo.PixelSpacing';
            end
            dcmdata(j).PixelSpacing = val;
            val = [];
            if isfield(tinfo,'SpacingBetweenSlices')
                val = tinfo.SpacingBetweenSlices;
            elseif isfield(tinfo,'SliceThickness')
                val = tinfo.SliceThickness;
            end
            dcmdata(j).SlcThk = val;
            val = '';
            if isfield(tinfo,'StudyDescription')
                val = tinfo.StudyDescription;
            end
            dcmdata(j).StudyDescription = val;
            if isfield(tinfo,'DiffusionBValue') && (length(tinfo.DiffusionBValue)==1)
                val = ['b',num2str(tinfo.DiffusionBValue)];
            elseif isfield(tinfo,'SeriesDescription') && ~isempty(tinfo.SeriesDescription)
                val = tinfo.SeriesDescription;
            elseif isfield(tinfo,'ProtocolName') && ~isempty(tinfo.ProtocolName)
                val = tinfo.ProtocolName;
            else val = tinfo.Modality;
            end
            dcmdata(j).Label = {val};
            val = '';
            if isfield(tinfo,'PatientID')
                val = tinfo.PatientID;
            end
            dcmdata(j).PatientID = val;
            dcmdata(j).img = [];
            if all(isfield(tinfo,{'Rows','Columns'}))
                dcmdata(j).d = [tinfo.Rows,tinfo.Columns];
            end
        end
        k = size(dcmdata(j).SlicePos,1)+1;
        if ~isfield(tinfo,'AcquisitionNumber') || isempty(tinfo.AcquisitionNumber)
            val = 0;
        else
            val = tinfo.AcquisitionNumber;
        end
        dcmdata(j).AcquisitionNumber(k) = val;
        val = nan(1,3);
        if isfield(tinfo,'ImagePositionPatient')
            val = tinfo.ImagePositionPatient;
        elseif isfield(tinfo,'SliceLocation')
            val(3) = tinfo.SliceLocation;
        else
            val(3) = ifn;
        end
        dcmdata(j).SlicePos(k,:) = val;
        val = [1,0,0,0,1,0]; % X-Y plane
        if isfield(tinfo,'ImageOrientationPatient')
            val = tinfo.ImageOrientationPatient;
        end
        dcmdata(j).SliceOrient(k,:) = val;

        % Read slice data from file:
        ySlope = 1;
        if isfield(tinfo,'RescaleSlope') && ~isempty(tinfo.RescaleSlope)
            ySlope = tinfo.RescaleSlope;
        elseif isfield(tinfo,'MRScaleSlope')
            ySlope = 1/tinfo.MRScaleSlope;
        end
        yInt = 0;
        if isfield(tinfo,'RescaleIntercept') && ~isempty(tinfo.RescaleIntercept)
            yInt = tinfo.RescaleIntercept;
        elseif isfield(tinfo,'MRScaleIntercept')
            yInt = -tinfo.MRScaleIntercept/ySlope;
        end
        timg = ySlope * double(dicomread(tinfo)) + yInt;
        if (length(tinfo.ImageType)>5) && strcmp(tinfo.ImageType(end-5:end),'MOSAIC')
            ns = int32(tinfo.LocationsInAcquisition(1));
            nslab = ceil(sqrt(double(ns)));
            md = size(timg);
            d = md/nslab;
            mosi = timg;
            timg = zeros([d,ns]);
            ci = reshape((1:nslab)'*ones(1,nslab),1,[]);
            ri = reshape(ones(nslab,1)*(1:nslab),1,[]);
            for i = 1:ns
                timg(:,:,i) = mosi( (ri(i)-1)*d(1)+(1:d(1)) , (ci(i)-1)*d(2)+(1:d(2)) );
            end
            dcmdata(j).d = d;
            kk = k - 1 + (1:ns);
            dcmdata(j).AcquisitionNumber(kk) = dcmdata(j).AcquisitionNumber(k);
            % Correct slice locations:
            RS = reshape(tinfo.ImageOrientationPatient,[3,2]);
            t = dcmdata(j).SlicePos(k,:)' + RS*(double(md)-double(d))'/2;
            RS(:,3) = cross(RS(:,1),RS(:,2));
            th = double(0:ns-1)*dcmdata(j).SlcThk;
            dcmdata(j).SlicePos(kk,:) = repmat(t',[ns,1]) + (RS*[zeros(ns,2),th']')';
            k = kk;
        end
        if isempty(dcmdata(j).img)
            dcmdata(j).img = timg;
        else
            dcmdata(j).img(:,:,k) = timg;
        end
    end
end
delete(hp);
catch err
    assignin('base','err',err);
    delete(hp);
    dcmdata = [];
    error('Error in readDICOM.m - Check variable "err" for more info.');
end

if ~isempty(dcmdata)
    
    % Separate by orientation:
    ndd = length(dcmdata);
    for i = 1:length(dcmdata)
        [SO,xi] = unique(dcmdata.SliceOrient,'rows');
        N = size(SO,1);
        if N~=1
            ii = [i,ndd+(1:N-1)];
            tdd = dcmdata(i);
            for j = 1:N
                ind = xi==i;
                dcmdata(ii(j)).img = tdd.img(:,:,:,ind);
                dcmdata(ii(j)).AcquisitionNumber = tdd.AcquisitionNumber(ind);
                dcmdata(ii(j)).SlicePos = tdd.SlicePos(ind,:);
                dcmdata(ii(j)).SliceOrient = SO(i,:);
            end
        end
        dcmdata(1).SliceOrient = SO(1,:);
    end

    % Sort by slice location
    for i = 1:length(dcmdata)
        snorm = cross(dcmdata(i).SliceOrient(1:3),dcmdata(i).SliceOrient(4:6));
        nslc = size(dcmdata(i).SlicePos,1);
        
        % Remove slices with no position:
        ind = any(isnan(dcmdata(i).SlicePos),2);
        if any(ind)
            dcmdata.img(:,:,ind) = [];
            dcmdata.SlicePos(ind,:) = [];
            dcmdata.AcquisitionNumber(ind) = [];
        end
        
        % Find repeated slice locations:
        [upos,pind] = unique(dcmdata(i).SlicePos,'rows');
        nupos = size(upos,1);
        if nslc==nupos
            ii = i;
        elseif mod(nslc,nupos)
            % Drop overlapping slices
            ii = i;
            rmslc = setxor(pind,1:nslc);
            dcmdata(i).img(:,:,rmslc) = [];
            dcmdata(i).AcquisitionNumber(rmslc) = [];
            dcmdata(i).SlicePos(rmslc,:) = [];
        else
            tdd = dcmdata(i);
            % Separate by AcquisitionNumber:
            a = tdd.AcquisitionNumber;
            ua = unique(a);
            nua = length(ua);
            if nua==nslc
                fprintf('Separating by file order\n')
                a = tdd.SlicePos;
                [~,ia] = unique(a,'rows');
                ii = [];
                while ~isempty(a)
                    if isempty(ii)
                        j = i;
                    else
                        j = length(dcmdata)+1;
                    end
                    ii = [ii,j];
                    dcmdata(j) = tdd;
                    dcmdata(j).img = tdd.img(:,:,:,ia);
                    dcmdata(j).AcquisitionNumber = tdd.AcquisitionNumber(ia);
                    dcmdata(j).SlicePos = tdd.SlicePos(ia,:);
                end
            else
                fprintf('Separating by AcquisitionNumber\n')
                ndd = length(dcmdata);
                ii = [i,(1:nua-1)+ndd];
                for j = 1:nua
                    ia = a==ua(j);
                    dcmdata(ii(j)) = tdd;
                    dcmdata(j).img = tdd.img(:,:,:,ia);
                    dcmdata(j).AcquisitionNumber = tdd.AcquisitionNumber(ia);
                    dcmdata(j).SlicePos = tdd.SlicePos(ia,:);
                end
            end
        end
        
        % Each slice location relative to first slice along normal vector:
        for j = ii;
            nslc = size(dcmdata(j).SlicePos,1);
            ds = dot(dcmdata(i).SlicePos - repmat(dcmdata(i).SlicePos(1,:),nslc,1),repmat(snorm,nslc,1),2);
            [ds,ids] = sort(ds);
            dcmdata(j).img = dcmdata(j).img(:,:,ids);
            dcmdata(j).AcquisitionNumber = dcmdata(j).AcquisitionNumber(ids);
            dcmdata(j).SlicePos = dcmdata(j).SlicePos(ids(1),:);
            dcmdata(j).PixelSpacing(3) = ds(2)-ds(1);
            
            % Check for varying slice spacing:
            dds = round(diff(ds),3);
            uds = unique(dds);
            if (numel(uds)>1) || (dds(1)>dcmdata(j).SlcThk)
                answer = questdlg('How would you like to compile these slices?',...
                    'Gapped acquisition','Concatenate','Insert Gaps','Cancel','Concatenate');
                switch answer
                    case 'Insert'
                        fval = str2double(inputdlg('Blank Slice Value:','',1,{num2str(min(dcmdata(j).img(:)))}));
                        ind = [1,cumsum(round(dds/min([uds,dcmdata(j).SlcThk])))+1];
                        timg = dcmdata(j).img;
                        dcmdata(j).img = ones([dcmdata(j).d,max(ind)])*fval;
                        dcmdata(j).img(:,:,ind) = timg;
                        dcmdata(j).PixelSpacing(3) = [];
                    case 'Cancel'
                         return;
                end
            end
            dcmdata(j).d(3) = size(dcmdata(j).img,3);
            dcmdata(j).PixelSpacing(3) = (ds(end)-ds(1))/double((dcmdata(j).d(3)-1));
        end
    end
    
    % User GUI to select single acquisition from set:
    ndd = length(dcmdata);
    if (ndd>1)
        if size(unique(reshape([dcmdata.d],[],3),'rows'),1)
            % Concatenate in 4D:
            dcmdata(1).img = cat(4,dcmdata(:).img);
            dcmdata(1).d(4) = size(dcmdata(1).img,4);
            dcmdata(1).Label = [dcmdata(:).Label];
            dcmdata(1).AcquisitionNumber = [dcmdata.AcquisitionNumber];
            dcmdata(1).SeriesNumber = {dcmdata.SeriesNumber};
            dcmdata(1).TemporalPositionIdentifier = [dcmdata.TemporalPositionIdentifier];
            dcmdata(1).DiffusionNumber = [dcmdata.DiffusionNumber];
            dcmdata(1).DiffusionDir = [dcmdata.DiffusionDir];
            dcmdata(1).TE = [dcmdata.TE];
            dcmdata(1).TR = [dcmdata.TR];
            dcmdata(2:end) = [];
        else
            str = strcat('(',cellfun(@(x)num2str(size(x,1)),{dcmdata(:).SlicePos},...
                                     'UniformOutput',false),...
                         ' Slices)',[dcmdata(:).Label]);
            % Show MIP preview for each:
            hf = figure('Colormap',gray);
            nrow = round(sqrt(ndd));
            ncol = ceil(ndd/nrow);
            for imont = 1:ndd
                subplot(nrow,ncol,imont),
                imshow(squeeze(max(dcmdata(imont).img,[],3)),[]);
                title(str{imont});
            end
            answer = listdlg('ListString',str,'SelectionMode','single',...
                             'ListSize',[300,300]);
            delete(hf);
            if isempty(answer)
                return % User cancelled the load
            else
                dcmdata = dcmdata(answer);
            end
        end
    end
    
    % Check for dimension match
    [d(1),d(2),d(3),d(4)] = size(dcmdata.img);
    if ~isempty(din) && ~all(d(1:3)==din)
        error(['Dimensions must match!',sprintf(' %f',d)])
    end

    % Prepare data to return:
    img = permute(dcmdata.img,[2,1,3,4]);
    dcmdata = rmfield(dcmdata,'img');
    [d(1),d(2),d(3),d(4)] = size(img);
    fov = dcmdata.PixelSpacing .* d(1:3);
    label = dcmdata.Label;
end
