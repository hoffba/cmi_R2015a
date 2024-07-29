function fn = brainReg(procdir,ref,hom,warpflag)
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
    fn_seg = fullfile(path_seg,[extractBefore(fn_seg,'.nii'),'.SynthSeg.nii.gz']);
    if isfile(fn_seg)
        ref.seg = readNIFTI(fn_seg);
    else
        ref.seg = ref.img > graythresh(ref.img);
    end
end
if ischar(hom)
    [img,label,fov,orient,info] = readNIFTI(hom);
    hom = struct('img',img,'label',label,'fov',fov,'orient',orient,'name',hom,'info',info);
end
if nargin<4
    warpflag = false;
end

% Initialize Reg Object
regObj = RegClass(false);
regObj.cmiObj(1).setImg(ref.img,ref.label,ref.fov,ref.orient,ref.name);
[~,homname] = fileparts(hom.name);
homname = regexprep(extractBefore(homname,'.nii'),'\.','_');
regObj.cmiObj(2).setImg(hom.img,hom.label,hom.fov,hom.orient,homname);
hom.img = [];

% Generate mask for registration (assumes everthing registered to HR-T2w)
mask = imsegfmm(ref.img,ref.img>prctile(ref.img,90,'all'),0.001);
ref.img = [];
regObj.cmiObj(1).img.mask.merge('replace',mask);
regObj.cmiObj(1).morphMask('d',[3 0]);
regObj.cmiObj(1).fillHoles;
regObj.cmiObj(1).morphMask('e',[3 0]);

% Determine output directory
odir = fullfile(procdir,['elxreg_',hom.label]);

% Align center of image as initial Translation transform
regObj.pts2M;

% Affine transform for first step
regObj.addElxStep('Affine','NumberOfResolutions',1,...
                           'FixedImagePyramidSchedule',ones(1,3),...
                           'MovingImagePyramidSchedule',ones(1,3),...
                           'MaximumNumberOfIterations',2000,...
                           'SP_A',50,...
                           'SP_a',200,...
                           'NumberOfSpatialSamples',5000,...
                           'DefaultPixelValue',0,...
                           'WriteResultImage','false');
if warpflag
    % Warping transform for longitudinal registrations
    fips = [max(round(8./regObj.cmiObj(1).img.voxsz),1),...
            max(round(3./regObj.cmiObj(1).img.voxsz),1)];
    mips = [max(round(8./regObj.cmiObj(2).img.voxsz),1),...
            max(round(3./regObj.cmiObj(2).img.voxsz),1)];
    gsp = [max(round(10./regObj.cmiObj(2).img.voxsz),1),...
           max(round( 3./regObj.cmiObj(2).img.voxsz),1)];
    % Warping registration between time points
    regObj.addElxStep('Warp','NumberOfResolutions',2,...
                         'FixedImagePyramidSchedule',fips,...
                         'MovingImagePyramidSchedule',mips,...
                         'MaximumNumberOfIterations',[2000,2000],...
                         'SP_A',50,...
                         'SP_a',[20000,2000],...
                         'NumberOfSpatialSamples',[5000,5000],...
                         'GridSpacingSchedule',gsp,...
                         'FinalGridSpacingInVoxels',gsp(4:6),...
                         'TransformBendingEnergy',10,...
                         'DefaultPixelValue',0,...
                         'WriteResultImage','false');
end

if ~isfolder(odir)
    mkdir(odir);
end
regObj.setOdir(odir);

regObj.setWait(true);
[stat,cmdstr] = regObj.runElx(false);

% Return name of registered image
fn = fullfile(odir,[homname,'_R.nii']);
