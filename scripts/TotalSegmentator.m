function seg = TotalSegmentator(ct,info,id,savepath)
% Whole body CT segmentation using TotalSegmentator


% First make sure TotalSegmentator is installed
TSpath = findTS(false);
if isempty(TSpath) || ~isfile(TSpath)
    return;
end

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
    % Will save a temporary image file for TotalSegmentator to read
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


function TSpath = findTS(flag)
    TSpath = [];
    if ispc
        str = 'where';
    else
        str = 'which';
    end
    [~,tpath] = system([str,' TotalSegmentator']);
    tpath = strsplit(tpath,'\n');
    ind = find(contains(tpath,'Python'),1,'last');
    if isempty(ind)
        if flag
            TSpath = 'TotalSegmentator not available. See documentation for installation instructions';
        else
            % Try to install TotalSegmentator
            stat = TSinstall;
            if stat
                TSpath = findTS(true);
            end
        end
    else
        TSpath = tpath{ind};
    end


function stat = TSinstall
    if ispc
        connect_str = ' && ';
    else
        connect_str = ' ; ';
    end
    cmd = {'pip install TotalSegmentator',...
           'pip install cupy-cuda11x cucim'};
    [stat,cmd_out] = system(strjoin(cmd,connect_str));
    




