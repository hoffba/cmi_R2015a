function stat = dcm2fld
stat = true;
savdir = '/mnt/cmi/projects/Skeletal_Studies/clin_bone/UMCC_2011_069_BoneMets/RawData/';

% Find and load all DICOM files in selected directory
[fname,fpath] = uigetfile('*','Select DICOM file:');
[~,~,ext] = fileparts(fname);
fnames = dir([fpath,filesep,'*',ext]);
fnames = {fnames(:).name};
% account for some dicom files with no extension
if isempty(ext)
    fnames2 = dir([path,filesep,'*.*']);
    fnames2 = {fnames2(:).name};
    fnames = setdiff(fnames,fnames2);
end
nf = length(fnames);

% Needed header fields:
tdata = struct;
tfields = {'Rows','Columns','PixelSpacing','Modality',...
    'TemporalPositionIdentifier','NumberOfFrames','RescaleSlope',...
    'RescaleIntercept','SliceLocation','InstanceNumber','ImagePositionPatient',...
    'AcquisitionNumber','OverlayData_0','SliceThickness','SpacingBetweenSlices'};
hp = waitbar(0,'','WindowStyle','modal',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)'); % Option to cancel load midway
assignin('base','hp',hp);
setappdata(hp,'canceling',0)
for i = 1:nf
% Check for Cancel button press
    if getappdata(hp,'canceling')
        return;
    end
    waitbar((i-1)/nf,hp,fnames{i});
    
    info = dicominfo(fullfile(fpath,fnames{i}));
    fstrs = tfields(isfield(info,tfields));
    for j = fstrs
        tdata(i).(j{1}) = info.(j{1});
    end
    tdata(i).mat = dicomread(fullfile(fpath,fnames{i}));
end
delete(hp)

% Determine data order:
if isfield(tdata,'OverlayData_0') % remove overlays
    tdata(~isempty([tdata(:).OverlayData_0])) = [];
end
nf = length(tdata);
acqn = ones(nf,1);
if isfield(tdata,'AcquisitionNumber')
    acqn = [tdata(:).AcquisitionNumber]';
    acqID = unique(acqn);
end
tpos = ones(nf,1);
if isfield(tdata,'TemporalPositionIdentifier')
    tpos = [tdata(:).TemporalPositionIdentifier]';
end
slcloc = ones(nf,1);
if isfield(tdata,'SliceLocation')
    slcloc = [tdata(:).SliceLocation]';
end
[~,ix] = sortrows([acqn,tpos,slcloc]);
acqn = acqn(ix); tpos = tpos(ix); slcloc = slcloc(ix);
tdata = tdata(ix);

% Save each acquisition separately
for i = 1:length(acqID)
    if isfield(tdata,'AcquisitionNumber')
        ind = (acqn==acqID(i));
    else
        ind = ones(1,nf);
    end
    ntot = sum(ind);
    ii = find(ind,1);
    
    % Check for time-series data
    tchk = true; t = 1; nt = 1;
    if isfield(tdata,'TemporalPositionIdentifier')
        t = unique(tpos(ind));
        nt = length(t);
        tchk = hist(tvec,t);
        tchk = all(tchk==tchk(1));
    end
    
    % Check 2D dimensions
    d1 = unique([tdata(ind).Rows]);
    d2 = unique([tdata(ind).Columns]);
    
    % Determine 3D slices
    ns = length(unique(slcloc(ind)));
    
    % Make sure matrix sizes are the same
    if (length(d1)==1) && (length(d2)==1) && tchk ...
            && (prod([ns,nt])==ntot)
        
        % Put slices in correct order
        img = reshape([tdata(ind).mat],[d1,d2,ns,nt]);
        
        % Calculate image FOV
        fov = [ double([d1,d2]).*tdata(ii).PixelSpacing(:)' , ...
                ns*tdata(ii).SliceThickness ];
            
        % Determine image labels
        if nt==1
            label = {tdata(ii).Modality};
        else
            label = cellfun(@num2str,num2cell(t),'UniformOutput',false);
        end

        % Now save the data as FLD
        fnout = fullfile(savdir,[info.PatientID,'_',info.StudyID,'_',...
            num2str(info.SeriesNumber),'_acq',num2str(acqID(i)),'.fld']);
        stat = saveFLD(fnout,img,label,fov,[]);
    end
end

