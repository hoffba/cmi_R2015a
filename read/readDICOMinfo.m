function info = readDICOMinfo(fpath)

info = struct('native_info',[],'label','','voxsz',[],'d',[],'fov',[],'orient',[]);

% Find DICOM files to load
fname = dir(fpath);
fname([fname.isdir]) = [];
fname = {fname(3:end).name};
ind = find(strcmpi('dicomdir',fname),1);
if ind
    % DICOM names are stored in "DICOMDIR" file
    %     hp = waitbar(0,'Reading DICOMDIR ...');
    fprintf('\nReading DICOMDIR ...\n');
    tinfo = dicominfo(fullfile(fpath,fname{ind}),'UseDictionaryVR',true);
    fname = fieldnames(tinfo.DirectoryRecordSequence);
    nf = length(fname);
    ind = false(1,nf);
    fprintf('Finding file locations ');
    for i = 1:nf
        if isfield(tinfo.DirectoryRecordSequence.(fname{i}),'ReferencedFileID')
            tdir = tinfo.DirectoryRecordSequence.(fname{i}).ReferencedFileID;
            fname{i} = strrep(tdir,'\',filesep);
        else
            ind(i) = true;
        end
        %         waitbar(i/nf,hp,'Finding file locations ...');
        if ~mod(i,round(nf/20))
            fprintf('.');
        end
    end
    %     delete(hp);
    fprintf('\n');
    fname(ind) = [];
else
    % All DICOMs are in single directory
    fname(~cellfun(@isdicom,fullfile(fpath,fname))) = [];

    % Files to exclude:
    exclude_fn = {'tdf'};
    fname(ismember(fname,exclude_fn)) = [];
end
nf = length(fname);

% Initialize structure of DICOM info
dcmdata = struct(...
                 'SeriesInstanceUID',[],...
                 'SeriesNumber',[],...
                 'd',[],...
                 'StudyDescription','',...
                 'PatientID','',...
                 'Label','', ...
                 'Orient',[],...
                 'TemporalPositionIdentifier',[],...
                 'DiffusionNumber',[],...
                 'DiffusionDir',[],...
                 'TE',[],...
                 'TR',[],...
                 'SlcThk',[],...
                 'PixelSpacing',[], ...
                 'SamplesPerPixel',[],...
             ... Lists for the following:
                 'fn',{},...
                 'AcquisitionNumber',[],...
                 'SlicePos',[]...
                 );
fn_ind = nan(nf,1); % Track which dcmdata index each file belongs to

% Loop over files
% * Images are sorted into 3D via SlcLoc,
%   and 4D via arrayed values:
%       MRI - TE, TR, ...
%       CT  - kV, ...
fprintf(['Loading DICOM info from %u files: ',repmat(' ',1,14)],nf);
for ifn = 1:nf
    fprintf([repmat('\b',1,14),'% 5u ( % 3u%% )'],ifn,round(ifn/nf*100));
    
    % Load DICOM file info:
    tinfo = dicominfo(fullfile(fpath,fname{ifn}),'UseDictionaryVR',true);
    if ifn==1 % Check manufacturer to read private tags:
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
        tinfo = dicominfo(fullfile(fpath,fname{ifn}),'UseDictionaryVR',true);
    end

    if isfield(tinfo,'Modality')

        [j,dcmdata] = find_index(dcmdata,tinfo);
        fn_ind(ifn) = j;

        % Add values to existing lists
        k = size(dcmdata(j).fn,1)+1;

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

    end
end
fprintf('\n');

n = length(dcmdata);
if n
    
    % Sort by slice location
    for i = 1:n
        [~,ix] = sort(dcmdata(i).SlicePos(:,3));
        dcmdata(i).fn = dcmdata(i).fn(ix);
        dcmdata(i).SlicePos = dcmdata(i).SlicePos(ix,:);
        dcmdata(i).AcquisitionNumber = dcmdata(i).AcquisitionNumber(ix);
    end
    
    % If multiple images per position, try and separate based on AcquisitionNumber
    ct = 0;
    for i = 1:n
        ii = i+ct;
        ns = size(dcmdata(ii).SlicePos,1);
        uS = unique(dcmdata(ii).SlicePos,'rows');
        nuS = length(uS);
        uA = unique(dcmdata(ii).AcquisitionNumber);
        nuA = length(uA);
        if (nuA>1) && (nuA<ns) && (nuS~=ns)
            if (nuA>1)
                disp('Separating images by AcquisitionNumber');
                A = dcmdata(ii).AcquisitionNumber;
                uA = unique(A);
                nuA = length(uA);
                dcmdata = [ dcmdata(1:(ii-1)),...
                    repmat(dcmdata(ii),1,nuA),...
                    dcmdata((ii+1):end) ];
                for j = 1:nuA
                    ind = ~(A==uA(j));
                    dcmdata(ii-1+j).fn(ind) = [];
                    dcmdata(ii-1+j).SlicePos(ind,:) = [];
                    dcmdata(ii-1+j).AcquisitionNumber(ind) = [];
                end
                ct = ct + nuA - 1;
            else
                disp('Separating images based on file order.');
                [~,ind] = sortrows(dcmdata(ii).SlicePos);
                nsep = ns/nuS;
                dcmdata = [ dcmdata(1:(ii-1)),...
                    repmat(dcmdata(ii),1,nsep),...
                    dcmdata((ii+1):end) ];
                for j = 1:nsep
                    jj = ii-1+j;
                    indj = ind(j:nsep:end);
                    dcmdata(jj).fn = dcmdata(jj).fn(indj);
                    dcmdata(jj).SlicePos = dcmdata(jj).SlicePos(indj,:);
                    dcmdata(jj).AcquisitionNumber = dcmdata(jj).AcquisitionNumber(indj);
                end
            end
        end
    end
    
    % Check if each slice is a separate acquisition
    uN = unique(cellfun(@(x)size(x,1),{dcmdata(:).SlicePos})); % Number of slices per set
    if (n>1) && (isscalar(uN)) && ismember(uN,1:2)
        % 1: Case for single-slice gapped CT data
        % 2: Case for 2-slice gapped CT data
        dcmdata(1).SlicePos = cat(1,dcmdata(:).SlicePos);
        dcmdata(1).AcquisitionNumber = cat(1,dcmdata(:).AcquisitionNumber);
        dcmdata(1).fn = cat(1,dcmdata(:).fn);
        dcmdata(2:end) = [];
    end

end

% Sort slices by location:
uind = dcmdata.SlicePos(:,3);
[~,ix,ic] = unique(uind,'rows');
ns = length(ix);
if (ns~=length(uind)) && isscalar(unique(histcounts(ic)))
    % MR Solutions DICOMs have no way to separate out b-values of DWI
    n4d =  d(3)/ns;
    uind = repmat(0:n4d-1,ns,1);
    uind = uind(:);
    ix = repmat(ix,n4d,1) + uind;

    dcmdata.SlicePos = dcmdata.SlicePos(ix,:);
    dcmdata.fn = [];
else
    % Takes the first slice image at a unique position
    dcmdata.fn = dcmdata.fn(ix);
    dcmdata.SlicePos = dcmdata.SlicePos(ix,:);
    dcmdata.AcquisitionNumber = dcmdata.AcquisitionNumber(:,ix);
end
if (n4d>1) && (d(3)==(n4d*nnz(uind==1)))
    ns = d(3)/n4d;
    dcmdata.img = reshape(dcmdata.img(:,:,ix),[d(1:2),ns,n4d]);
    dcmdata.SlicePos = dcmdata.SlicePos(ix(1:ns),:);
    % Determine labeling:
    dcmdata.Label = cell(1,n4d);
    ind = ix(1:ns:d(3));
    if isfield(dcmdata,'TR') && length(unique(dcmdata.TR))>1
        dcmdata.Label = strcat(dcmdata.Label,'TR=',...
            cellfun(@num2str,num2cell(dcmdata.TR(ind)),...
            'UniformOutput',false),'; ');
    end
    if isfield(dcmdata,'TE') && length(unique(dcmdata.TE))>1
        dcmdata.Label = strcat(dcmdata.Label,'TE=',...
            cellfun(@num2str,num2cell(dcmdata.TE(ind)),...
            'UniformOutput',false),'; ');
    end
    ind = cellfun(@isempty,dcmdata.Label);
    if any(ind)
        dcmdata.Label(ind) = cellfun(@num2str,num2cell(find(ind)),'UniformOutput',false);
    end
end

% Calculate 3D orientation matrix (compatible with NIfTI orientation)
info.orient = [ [  dcmdata.Orient(1:3)*dcmdata.PixelSpacing(1) ,...
                   dcmdata.Orient(4:6)*dcmdata.PixelSpacing(2) ,...
                   (dcmdata.SlicePos(end,:)-dcmdata.SlicePos(1,:))'/(d(3)-1) ,...
                   dcmdata.SlicePos(1,:)' ] ;...
                   0 0 0 1];


% Prepare data to return:
info.native_info.format = 'DICOM';

info.fov = [dcmdata.PixelSpacing',1] .* d(1:3);
if size(dcmdata.SlicePos,1)>1
    info.fov(3) = abs(sqrt(sum(diff(dcmdata.SlicePos([1,end],:),1).^2,2))) * d(3)/(d(3)-1);
elseif ~isempty(dcmdata.SlcThk)
    info.fov(3) = d(3) * dcmdata.SlcThk;
end

info.label = dcmdata.Label;



function [j,dcmdata] = find_index(dcmdata,tinfo)

    % Check if Series already exists in structure:
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
    elseif isfield(tinfo,'B_value') && (isscalar(tinfo.B_value))
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
        if ~isempty(tinfo.DiffusionDirectionX) && (isscalar(tinfo.DiffusionDirectionX))
            diffD(1) = tinfo.DiffusionDirectionX;
        end
        if ~isempty(tinfo.DiffusionDirectionY) && (isscalar(tinfo.DiffusionDirectionY))
            diffD(2) = tinfo.DiffusionDirectionY;
        end
        if ~isempty(tinfo.DiffusionDirectionZ) && (isscalar(tinfo.DiffusionDirectionZ))
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
    d = [tinfo.Rows,tinfo.Columns];

    j = [];
    if ~isempty(dcmdata)
        j = find( strcmp(tinfo.SeriesInstanceUID,{dcmdata(:).SeriesInstanceUID}) ...
            & (serN==[dcmdata(:).SeriesNumber]) ...
            & (tempN == [dcmdata(:).TemporalPositionIdentifier]) ...
            & (diffN == [dcmdata(:).DiffusionNumber]) ...
            & ismember([dcmdata(:).DiffusionDir]',diffD','rows')' ...
            & (TE==[dcmdata(:).TE]) ...
            & (TR==[dcmdata(:).TR]) ...
            & cellfun(@(x)all(x==[tinfo.Rows,tinfo.Columns]),{dcmdata.d}) ,1);
    end
    if isempty(j) % Initialize new series in structure:
        
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
            val = tinfo.PixelSpacing;
        end
        dcmdata(j).PixelSpacing = val;
        
        val = [];
        if isfield(tinfo,'SliceThickness')
            val = tinfo.SliceThickness;
        elseif isfield(tinfo,'SpacingBetweenSlices')
            val = abs(tinfo.SpacingBetweenSlices);
        end
        dcmdata(j).SlcThk = val;
        
        val = [1 0 0 0 1 0];
        if isfield(tinfo, 'ImageOrientationPatient')
            val = tinfo.ImageOrientationPatient;
        end
        dcmdata(j).Orient = val;
        
        val = '';
        if isfield(tinfo,'StudyDescription')
            val = tinfo.StudyDescription;
        end
        dcmdata(j).StudyDescription = val;
        
        if isfield(tinfo,'DiffusionBValue') && (isscalar(tinfo.DiffusionBValue))
            val = ['b',num2str(tinfo.DiffusionBValue)];
        elseif isfield(tinfo,'SeriesDescription') && ~isempty(tinfo.SeriesDescription)
            val = tinfo.SeriesDescription;
        elseif isfield(tinfo,'ProtocolName') && ~isempty(tinfo.ProtocolName)
            val = tinfo.ProtocolName;
        else
            val = tinfo.Modality;
        end
        dcmdata(j).Label = {val};
        
        val = '';
        if isfield(tinfo,'PatientIdentityRemoved') && strcmpi(tinfo.PatientIdentityRemoved,'yes') ...
                && isfield(tinfo,'PatientName')
            if isfield(tinfo.PatientName,'FamilyName')
                val = tinfo.PatientName.FamilyName;
            end
        elseif isfield(tinfo,'PatientID')
            val = tinfo.PatientID;
        elseif isfield(tinfo,'StudyID')
            val = tinfo.StudyID; % added this b/c some cases did not have PatientName tag in dicom
        end
        dcmdata(j).PatientID = val;
        
        val = '';
        if isfield(tinfo,'StudyDate')
            val = tinfo.StudyDate;
        end
        dcmdata(j).StudyDate = val;
        
        if all(isfield(tinfo,{'Rows','Columns'}))
            dcmdata(j).d = [tinfo.Rows,tinfo.Columns];
        end
    end