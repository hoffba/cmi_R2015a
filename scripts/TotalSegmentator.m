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
if ischar(ct)
    if isfile(ct) && contains(ct,'.nii')
        ct_fname = ct;
        savepath = fileparts(ct);
    else
        seg = 'Invalid input file name';
        return;
    end
elseif nargin==4
    cleanup_chk = true;
    ct_fname = fullfile(savepath,[id,'.nii']);
    saveNIFTI(ct_fname,ct,id,info.fov,info.orient)
else
    seg = 'Invalid inputs';
    return;
end

% Run TotalSegmentator
system(['cd /D ',savepath,' & ',...
        'python ',TSpath,...
        ' -i ',ct_fname,...
        ' -o ',savepath,...
        ' --ml']);
% Clean up
if cleanup_chk
    delete(ct_fname);
end
% Compile resulting lobe segmentations
seg = readNIFTI(fullfile(savepath,['TotalSegmentator_',id,'.nii']));
seg(~ismember(seg,13:17)) = 0; % Remove non-lung

