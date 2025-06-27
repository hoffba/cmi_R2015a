function Mreg = transformImage(fn_tf,flag_nn, M,voxsz,orient)

procdir = fileparts(fn_tf);

if ischar(M) && isfile(M)
    fn_flag = true;
    fn_M = M;
elseif isnumeric(M) && nargin==5
    % Write image to file
    fn_flag = false;
    d = size(M);
    fn_M = fullfile(procdir,'temp_image.nii');
    saveNIFTI(fn_M,M,{'temp'},voxsz.*d,orient);
else
    error('Invalid inputs')
end

p = readTransformParam(fn_tf);
p.ResultImageFormat = 'nii';
% Write new TransformParameter files if NN is needed
if flag_nn
    p.ResampleInterpolator = 'FinalNearestNeighborInterpolator';
    fn_tf = insertBefore(fn_tf,'.txt','.NN');
end
writeElxStruct2Txt(p,fn_tf);

% Transformix command
cmdstr = sprintf('transformix -in "%s" -out "%s" -tp "%s"',fn_M,procdir,fn_tf);
system(cmdstr);

if fn_flag
    % Zip result file
    fn_result = fullfile(procdir,'result.nii');
    gzip(fn_result);
    delete(fn_result);
    fn_result = [fn_result,'.gz'];
    % Rename result file
    fn = flip(extractBefore(flip(fn_M),filesep));
    fn = insertBefore(fn,'.nii','.reg');
    Mreg = fullfile(procdir,fn);
    movefile(fn_result,Mreg);
else
    % Return transformed image
    Mreg = readNIFTI(fullfile(procdir,'result.nii'));
    % Delete temp image file
    delete(fn_M);
end
