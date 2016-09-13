function [img,info] = readDICOM(fname,D)

img = []; info = [];
if nargin==1
    D = [];
end
if isdir(fname)
    fpath = fname;
    ext = '*';
else
    [fpath,fname,ext] = fileparts(fname);
end
if ~any(strcmp(ext,{'.dcm','.1'}))
    ext = [];
end

% Find and load all DICOM files in selected directory
if strcmpi(fname,'dicomdir')
    % DICOM names are stored in "DICOMDIR" file
    hp = waitbar(0,'Reading DICOMDIR ...');
    tinfo = dicominfo(fullfile(fpath,fname));
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
    nf = length(fnames);
end

if nf==0, return; end

dcmdata = struct('info',{},...
                 'SeriesInstanceUID',{},... % for sorting multiple images
                 'img',{},...
                 'd',{},...
                 'AcquisitionNumber',{},...
                 'SeriesNumber',{},...
                 'TemporalPositionIdentifier',{},...
                 'DiffusionBValue',{},...
                 'RepetitionTime',{},...
                 'EchoTime',{},...
                 'SlicePos',{},...
                 'SlcThk',{},...
                 'PixelSpacing',{},...
                 'StudyDescription',{},...
                 'PatientID',{},...
                 'Label',{},...
                 'Time',{});
             
% Initialize sorting options (in addition to slice position):
sortC = {'SeriesInstanceUID',           0;...
         'AcquisitionNumber',           1;...
         'SeriesNumber',                1;...
         'TemporalPositionIdentifier',  1;...
         'DiffusionBValue',             1;...
         'RepetitionTime',              1;...
         'EchoTime',                    1};
sortN = size(sortC,1);

% Read in all DICOM slices:
hp = waitbar(0,'','WindowStyle','modal',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)'); % Option to cancel load midway
setappdata(hp,'canceling',0);
ndcm = 0;
for ifn = 1:nf
    
    % Check for Cancel button press
    if getappdata(hp,'canceling')
        img = []; info = [];
        delete(hp);
        return;
    end
    waitbar((ifn-1)/nf,hp,fnames{ifn});
    
    % Load DICOM file info:
    tinfo = dicominfo(fullfile(fpath,fnames{ifn}),'UseDictionaryVR',true);
    
    % Check dimensions against desired dimensions:
    nD = [tinfo.Rows,tinfo.Columns];
    if isempty(D) || all(nD==D(1:2))
    
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
            tinfo = dicominfo(fullfile(fpath,fnames{ifn}),'UseDictionaryVR',true);
        end

        % Make sure it's an image with desired dimensions before proceeding
        if isfield(tinfo,'Modality') && (isempty(D) ...
                || ((tinfo.Rows==D(1)) && (tinfo.Columns==D(2))))

            % Check if Series already exists in structure:
            ichk = zeros(ndcm,sortN);
            for i = 1:sortN
                fld = sortC{i,1};
                if isfield(tinfo,fld) && ~isempty(tinfo.(fld))
                    if sortC{i,2} % Numerical
                        ichk(:,i) = (tinfo.(fld)==[dcmdata(:).(fld)]);
                    else          % String
                        ichk(:,i) = strcmp(tinfo.(fld),{dcmdata(:).(fld)});
                    end
                else
                    if sortC{i,2}
                        ichk(:,i) = isnan([dcmdata(:).(fld)]);
                    else
                        ichk(:,i) = strcmp('',{dcmdata(:).(fld)});
                    end
                end
            end
            j = find(all(ichk,2));

            if isempty(j)
                % Initialize new series in structure:
                ndcm = ndcm+1;
                j = ndcm;
                dcmdata(j).info = tinfo;
                for i = 1:sortN
                    fld = sortC{i,1};
                    if isfield(tinfo,fld) && ~isempty(tinfo.(fld))
                        val = tinfo.(fld);
                    elseif sortC{i,2}
                        val = nan;
                    else
                        val = '';
                    end
                    dcmdata(j).(fld) = val;
                end
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
                if isfield(tinfo,'DiffusionBValue') && ~isempty(tinfo.DiffusionBValue)
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
                dcmdata(j).d = nD;
            end
            k = size(dcmdata(j).SlicePos,1)+1;
            val = nan(1,3);
            if isfield(tinfo,'ImagePositionPatient')
                val = tinfo.ImagePositionPatient;
            elseif isfield(tinfo,'SliceLocation')
                val(3) = tinfo.SliceLocation;
            else
                val(3) = ifn;
            end
            dcmdata(j).SlicePos(k,:) = val;

            val = nan;
            if isfield(tinfo,'TimeAfterStart')
                val = tinfo.TimeAfterStart;
            end
            dcmdata(j).Time(k) = val;

            % Read slice data from file:
            if isfield(tinfo,'RescaleSlope')
                ySlope = tinfo.RescaleSlope;
            elseif isfield(tinfo,'MRScaleSlope')
                ySlope = 1/tinfo.MRScaleSlope;
            else ySlope = 1;
            end
            if isfield(tinfo,'RescaleIntercept')
                yInt = tinfo.RescaleIntercept;
            elseif isfield(tinfo,'MRScaleIntercept')
                yInt = -tinfo.MRScaleIntercept/ySlope;
            else yInt = 0;
            end
            timg = ySlope * double(dicomread(tinfo)) + yInt;
            dcmdata(j).img(:,:,k) = timg;
        end
    end
end
delete(hp);

if ndcm==0
    return;
end

% Determine labels for series data base on what's different:
lstr = repmat({''},ndcm,1);
tval = [dcmdata(:).AcquisitionNumber]';
if length(unique(tval))>1
    lstr = strcat(lstr,'a',cellfun(@num2str,num2cell(tval),'UniformOutput',false));
end
tval = [dcmdata(:).SeriesNumber]';
if length(unique(tval))>1
    lstr = strcat(lstr,'s',cellfun(@num2str,num2cell(tval),'UniformOutput',false));
end
tval = [dcmdata(:).TemporalPositionIdentifier]';
if length(unique(tval))>1
    lstr = strcat(lstr,'t',cellfun(@num2str,num2cell(tval),'UniformOutput',false));
end
tval = [dcmdata(:).RepetitionTime]';
if length(unique(tval))>1
    lstr = strcat(lstr,'tr',cellfun(@num2str,num2cell(tval),'UniformOutput',false));
end
tval = [dcmdata(:).EchoTime]';
if length(unique(tval))>1
    lstr = strcat(lstr,'te',cellfun(@num2str,num2cell(tval),'UniformOutput',false));
end

% Sort slices by SlicePosition
for i = 1:length(dcmdata)
    [~,ia,ic] = unique(fliplr(dcmdata(i).SlicePos),'rows');
    nn = cellfun(@(x)nnz(x==ic),num2cell(ia));
    if all(nn==1) % no position repeats
        dcmdata(i).SlicePos = dcmdata(i).SlicePos(ic,:);
        dcmdata(i).Time = dcmdata(i).Time(ic);
        dcmdata(i).img = dcmdata(i).img(:,:,ic);
    elseif all(nn==nn(1)) % can be split into 4D
        % Assumes slices were read in order
        dcmdata(i).SlicePos = dcmdata(i).SlicePos(ia,:);
        [~,ix] = sort(ic);
        dcmdata(i).Time = reshape(dcmdata(i).Time(ix),[],nn(1));
        dcmdata(i).img = reshape(dcmdata(i).img(:,:,ix),...
            dcmdata(i).d(1),dcmdata(i).d(2),[],nn(1));
        dcmdata(i).Label = strcat(dcmdata(i).Label,'_',...
            cellfun(@num2str,num2cell(1:nn(1))));
    else
        error('Cannot reconcile slices.');
    end
    dcmdata(i).info.ImagePositionPatient = dcmdata(i).SlicePos(1,:);
    dcmdata(i).Label = strcat(dcmdata(i).Label,'_',lstr{i});
end

gflag = 0; % Flag for number of slices per location group
oimg = [];
uN = unique(cellfun(@(x)size(x,1),{dcmdata(:).SlicePos}));
if (ndcm==1) && (length(unique(round(diff(dcmdata.SlicePos(:,3)),4)))==2)
% Case for 2-slice gapped CT data
    gflag = 2;
    oimg = dcmdata.img;
    oloc = dcmdata.SlicePos;
elseif (length(uN)==1) && ismember(uN,1:2)
% 1: Case for single-slice gapped CT data
% 2: Case for 2-slice gapped CT data
    gflag = size(dcmdata(1).SlicePos,1);
    oimg = cat(3,dcmdata(:).img);
    oloc = cat(1,dcmdata(:).SlicePos);
end
if gflag
    answer = questdlg('How would you like to compile these slices?',...
        'Gapped CT','Concatenate','Insert Gaps','Cancel','Concatenate');
    if ~strcmp(answer,'Cancel')
        tdata = dcmdata(1);
        if strcmp(answer,'Insert Gaps')
            fval = str2double(inputdlg('Blank Slice Value:','',1,{'-1024'}));
        % Determine gaps:
            [~,ix] = sort(oloc(:,3));
            dxyz = sqrt(sum(diff(oloc(ix,:),1).^2,2));
            d = size(oimg);
            if gflag == 2
                dnew = floor((dxyz(1)+dxyz(2))/dxyz(1));
                dz = abs(dxyz(1)+dxyz(2))/dnew;
                d(3) = d(3)*dnew/2;
                tmat = ones(d)*fval;
                ind = round(dnew/2):dnew:d(3);
                ind = [ind;ind+1];
            else % Single-slice
                dnew = floor(dxyz(1)/tdata.SlcThk);
                dz = abs(dxyz(1))/dnew;
                d(3) = d(3)*dnew;
                tmat = ones(d)*fval;
                ind = round(dnew/2):dnew:d(3);
            end
            disp(['Image slices are now: ',num2str(ind(:)')])
            tmat(:,:,ind(:)) = oimg(:,:,ix);
            tdata.SlicePos = (0:d(3)-1)'*dz + oloc(1,:);
            tdata.SlcThk = dz;
            tdata.img = tmat;
        else
            tdata.img = oimg;
            tdata.SlicePos = oloc;
        end
        dcmdata = tdata;
    end
end

% User GUI to select single acquisition from set:
if (length(dcmdata)>1)
    dd = cellfun(@size,{dcmdata(:).img},'UniformOutput',false);
    ld = cellfun(@length,dd);
    if all(ld==ld(1)) && all(cellfun(@(x)all(x==dd{1}),dd))
        dcmdata(1).img = cat(4,dcmdata(:).img);
        dcmdata(1).Label = [dcmdata(:).Label];
        dcmdata(1).Time = cat(2,dcmdata(:).Time);
        dcmdata(2:end) = [];
    else
        str = strcat('(',cellfun(@(x)num2str(size(x,1)),{dcmdata(:).SlicePos},...
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
end

ind = isnan(dcmdata.SlicePos(:,3));
dcmdata.img(:,:,ind) = [];
dcmdata.SlicePos(ind,:) = [];

% Check for 3D dimension match
[d(1),d(2),d(3),d(4)] = size(dcmdata.img);
if ~isempty(D) && ~all(d(1:3)==D)
    error(['Dimensions must match!',sprintf(' %f',d)])
end

% Prepare data to return:
img = dcmdata.img;
info = dcmdata.info;
info.Label = dcmdata.Label;
info.NumberOfSlices = size(img,3);
info.ArrayDim = size(img,4);
info.Time = dcmdata.Time;
info.FileName = fpath;
info.Format = 'DICOM';

