function seg = TotalSegmentator(ct,info,id,savepath)
% Whole body CT segmentation using TotalSegmentator

seg = [];

% First make sure TotalSegmentator is installed
if ispc
    str = 'where';
else
    str = 'which';
end
[~,TSpath] = system([str,' TotalSegmentator']);
TSpath = strsplit(TSpath,'\n');
ind = find(contains(TSpath,'Python'),1,'last');
if isempty(ind)
    seg = 'TotalSegmentator not found';
    return;
end
TSpath = TSpath{ind};

% Manage inputs
cleanup_chk = false;
if ischar(ct)
    if ~(isfile(ct) && contains(ct,'.nii'))
        seg = 'Invalid input file name';
        return;
    end
    ct_fname = ct;
    [savepath,id] = fileparts(ct);
    id = extractBefore(id,'.');
elseif nargin==4
    cleanup_chk = true;
    ct_fname = fullfile(savepath,[id,'.nii']);
    saveNIFTI(ct_fname,ct,id,info.fov,info.orient)
else
    seg = 'Invalid inputs';
    return;
end

% Run TotalSegmentator
sv_name = fullfile(savepath,[id,'.TotalSegmentator.nii.gz']);
system(['cd /D ',savepath,' & ',...
        'python ',TSpath,...
        ' -i ',ct_fname,...
        ' -o ',sv_name,...
        ' --ml']);
% Clean up
if cleanup_chk
    delete(ct_fname);
end
% Compile resulting lobe segmentations
seg = readNIFTI(sv_name);
seg(~ismember(seg,13:17)) = 0; % Remove non-lung

function stat = TSinstall
    'pip install TotalCommander'
    'pip install cupy-cuda11x cucim'