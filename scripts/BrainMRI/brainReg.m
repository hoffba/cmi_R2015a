function fn = brainReg(procdir,ref,hom,flag)
% Inputs:
%   procdir = Folder containing study data
%   ref = Reference image (Ax MPR):
%           struct created in dipg
%       OR  filename
%   hom = Homologous image to register to ref
%           struct created in dipg
%       OR  filname
%   flag = T: register between timepoints with warp
%          F: affine registration within study

if ischar(ref)
    [img,label,fov,orient,info] = readNIFTI(ref);
    ref = struct('img',img,'label',label,'fov',fov,'orient',orient,'name',ref,'info',info);
    % Find segmentation to load
    [path_seg,fn_seg] = fileparts(ref.name);
    fn_seg = fullfile(path_seg,[extractBefore(fn_seg,'.'),'.SynthSeg.nii.gz']);
    if isfile(fn_seg)
        ref.seg = readNIFTI(fn_seg);
    else
        error('Could not find SynthSeg')
    end
end
if ischar(hom)
    [img,label,fov,orient,info] = readNIFTI(hom);
    hom = struct('img',img,'label',label,'fov',fov,'orient',orient,'name',hom,'info',info);
end
if nargin<4
    flag = false;
end

regObj = RegClass(false);
regObj.cmiObj(1).setImg(ref.img,ref.label,ref.fov,ref.orient,ref.name);
regObj.cmiObj(1).img.mask.merge('replace',logical(ref.seg));
% Dilate mask by 5mm
regObj.UDpreproc('dilaten',[round(5./regObj.cmiObj(1).img.voxsz);0,0,0])
regObj.addElxStep('Affine','NumberOfResolutions',1,...
                           'FixedImagePyramidSchedule',ones(1,3),...
                           'MovingImagePyramidSchedule',ones(1,3),...
                           'MaximumNumberOfIterations',4000,...
                           'SP_A',50,...
                           'SP_a',2000,...
                           'NumberOfSpatialSamples',5000,...
                           'DefaultPixelValue',0,...
                           'WriteResultImage','false');
if flag
    % Warping registration between time points
    regObj.addElxStep('Warp','NumberOfResolutions',1,...
                         'FixedImagePyramidSchedule',ones(1,3),...
                         'MovingImagePyramidSchedule',ones(1,3),...
                         'MaximumNumberOfIterations',2000,...
                         'SP_A',50,...
                         'SP_a',10000,...
                         'NumberOfSpatialSamples',10000,...
                         'GridSpacingSchedule',5*ones(1,3),...
                         'TransformBendingEnergy',50,...
                         'DefaultPixelValue',0,...
                         'WriteResultImage','false');

end
[~,homname] = fileparts(hom.name);
homname = regexprep(extractBefore(homname,'.nii'),'\.','_');
regObj.cmiObj(2).setImg(hom.img,hom.label,hom.fov,hom.orient,homname);

% Set output directory
odir = fullfile(procdir,['elxreg_',hom.label]);
if ~isfolder(odir)
    mkdir(odir);
end
regObj.setOdir(odir);

regObj.setWait(true);
[stat,cmdstr] = regObj.runElx(false);

% Return name of registered image
fn = fullfile(odir,[homname,'_R.nii']);
