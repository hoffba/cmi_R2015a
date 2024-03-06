function seg = TotalSegmentator(ct,info,id,savepath)
% Whole body CT segmentation using TotalSegmentator
seg = [];

% Manage inputs
cleanup_chk = false;
if ischar(ct)
    if ~(isfile(ct) && contains(ct,'.nii'))
        seg = 'Invalid input file name';
        return;
    end
    [savepath,id] = fileparts(ct);
    id = extractBefore(id,'.');
    ct_fname = ct;
elseif nargin==4
    % Will save a temporary image file for TotalSegmentator to read
    cleanup_chk = true;
    ct_fname = fullfile(savepath,[id,'.TotalSegTEMP.nii']);
    saveNIFTI(ct_fname,ct,id,info.fov,info.orient)
else
    fprintf('Invalid inputs\n');
    return;
end

cmd = '';
scriptpath = fullfile(fileparts(which('cmi')),'shellscripts');
segname = fullfile(savepath,[id,'.TotalSegmentator.nii.gz']);
if ispc
    sh_path = fullfile(scriptpath,'run_TotalSegmentator_Win.cmd');
    cmd = ['"',sh_path,'" "',ct_fname,'" "',segname,'"'];
elseif isunix
    sh_path = fullfile(scriptpath,'run_TotalSegmentator_Linux.sh');
    % cmd = ['module load python3.10-anaconda/2023.03 ...' ...
    %       '& conda init bash ...' ...
    %       '& bash -i ' sh_path ' ...' ...
    %             ct_fname,' ',...
    %             segname];
    cmd = ['bash -i ' sh_path ' ' ...
                      ct_fname,' ',...
                      segname];
else
end

% Run TotalSegmentator
if ~isempty(cmd)
    system(cmd);
end

% Clean up
if cleanup_chk
    delete(ct_fname);
end
% Compile resulting lobe segmentations
if nargout
    seg = readNIFTI(segname);
    seg(~ismember(seg,13:17)) = 0; % Remove non-lung
end

    



