function j = lungreg_BH(exp,exp_info,exp_mask,ins,ins_info,ins_mask,elxdir,ID,qcheck,flag) 
%OLD --lungreg_BH(ID,elxdir,regObj,qcheck)

% ***_info = struct with fields:
%                         .name
%                         .label
%                         .fov
%                         .orient
%                         .d
%                         .voxsz
%                         .voxvol


% This function runs key steps for image RegObj0istration using Ben's RegObj0
% program. This program will execute a RegObj0istration using preloaded data,
% segmentation masks and settings. The automated features include: generate
% masks, perform match VOI, delete moving mask, separate lungs if available and execute queue.

% To execute follow these steps:
% 1. >> RegObj0
% 2. Set Options/Image Limits [-Inf 0] for both Homologous and Reference
% 3. Set Options/VOI Dilate [20 20 10]
% 4. Load Elastix/Elastix Parameters
% Steps 1-4 only have to be once.
% The following steps are repeated for each RegObj0istration case
% 1. Load fixed Image (RED)
% 2. Append fixed segmentation map
% 3. Move 4D slider to 1st position (Image)
% 4. Load moving Image (Blue)
% 5. Append moving segmentation map
% 6. Move 4D slider to 1st postion (Image)
% 7. >> batch_lungreg_v3(C,RegObj0,N)

% Input runs out of batch_process.m
% Fixed_dFile = file location of fixed data
% Moving_dFile = file location of moving data
% Fixed_mask = mask for fixed data
% Moving_mask = mask for moving data
% RegObj0 = reg object

if nargin<4
    qcheck = true;
else
    qcheck = logical(qcheck);
end
if nargin<10
    flag = false;
end

%% Set up Reg object:
regObj = RegClass(false);
regObj.cmiObj(1).setImg(exp,exp_info.label,exp_info.fov,exp_info.orient,exp_info.name);
regObj.cmiObj(1).img.mask.merge('replace',exp_mask);
regObj.cmiObj(2).setImg(ins,ins_info.label,ins_info.fov,ins_info.orient,ins_info.name);
regObj.cmiObj(2).img.mask.merge('replace',ins_mask);

%% Is the fixed image incremental or continuous?
nslc = regObj.cmiObj(1).img.dims(3);
gapchk = nnz(std(reshape(regObj.cmiObj(1).img.mat(:,:,:,1),[],nslc),[],1)) < nslc;
if gapchk
    disp('Fixed image is incremental')
end
nslc = regObj.cmiObj(2).img.dims(3);
gapchk(2) = nnz(std(reshape(regObj.cmiObj(2).img.mat(:,:,:,1),[],nslc),[],1)) < nslc;
if gapchk(2)
    disp('Moving image is incremental')
end

%% Generate Registration Masks
% mask = ismember(regObj.cmiObj(1).img.mat(:,:,:,2),1:60); % create fixed mask
% regObj.cmiObj(1).img.mask.merge('replace',mask); % set cmiObj mask
% 
% mask = ismember(regObj.cmiObj(2).img.mat(:,:,:,2),1:60);   % Create moving mask for Initial transform
% regObj.cmiObj(2).img.mask.merge('replace',mask);

%% Set up preprocessing options
regObj.UDpreproc('filtn',[3,3,0;3,3,0],'dilaten',[10,10,10;0,0,0],'clamp',[-1000,0;-1000,0])

%% Set up registration resolutions
f_sched = {8*ones(1,3), 10*ones(1,3), [5*ones(1,3),2*ones(1,3)]};
if gapchk(1)
    f_sched{1}(3) = 1;
    f_sched{2}(3) = 1;
    f_sched{3}([3,6]) = 1;
end
m_sched = {8*ones(1,3), 10*ones(1,3), [5*ones(1,3),2*ones(1,3)]};
if gapchk(2)
    m_sched{1}(3) = 1;
    m_sched{2}(3) = 1;
    m_sched{3}([3,6]) = 1;
end

%% Single res third step for faster reg:
nres = 2;
maxiter = [10000 2000];
SP_A = [50 50];
SP_a = [100000 25000];
samp = [5000 5000];
gridsp = [5*ones(1,3),2*ones(1,3)];
if flag
    nres = 1;
end

%% Add steps to the schedule
if ~isfolder(elxdir)
    mkdir(elxdir)
end
regObj.setOdir(elxdir);
regObj.addElxStep('Affine','NumberOfResolutions',1,...
                           'FixedImagePyramidSchedule',f_sched{1},...
                           'MovingImagePyramidSchedule',f_sched{1},...
                           'MaximumNumberOfIterations',4000,...
                           'SP_A',50,...
                           'SP_a',2000,...
                           'NumberOfSpatialSamples',2000,...
                           'ResultImageFormat','nii',...
                           'DefaultPixelValue',-2000);
                      
regObj.addElxStep('Warp','NumberOfResolutions',1,...
                         'FixedImagePyramidSchedule',f_sched{2},...
                         'MovingImagePyramidSchedule',m_sched{2},...
                         'MaximumNumberOfIterations',20000,...
                         'SP_A',50,...
                         'SP_a',200000,...
                         'NumberOfSpatialSamples',2000,...
                         'GridSpacingSchedule',10*ones(1,3),...
                         'TransformBendingEnergy',50,...
                         'ResultImageFormat','nii',...
                         'DefaultPixelValue',-2000);
                    
regObj.addElxStep('Warp','NumberOfResolutions',nres,...
                         'FixedImagePyramidSchedule',f_sched{3}(1:nres),...
                         'MovingImagePyramidSchedule',m_sched{3}(1:nres),...
                         'MaximumNumberOfIterations',maxiter(1:nres),...
                         'SP_A',SP_A(1:nres),...
                         'SP_a',SP_a(1:nres),...
                         'NumberOfSpatialSamples',samp(1:nres),...
                         'GridSpacingSchedule',gridsp(1:nres*3),...
                         'TransformBendingEnergy',25,...
                         'ResultImageFormat','nii',...
                         'DefaultPixelValue',-2000);

%% Run Reg Function then queue
disp([ID,': Initial transform:'])
regObj.VOI2Tform;
pause(1);

disp([ID,': Remove Moving Mask'])
regObj.cmiObj(2).clearMask;

disp([ID,': Enqueue Registration:'])
regObj.cmiObj(1).setVec(1);
regObj.cmiObj(2).setVec(1);
regObj.setWait(~qcheck);
regObj.runElx(qcheck);
pause(1);

delete(regObj);



