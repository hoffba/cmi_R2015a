function lungreg_BH(Fixed_dFile, Moving_dFile, Fixed_mFile, Moving_mFile, regObj)
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


%% Elxreg directory for results:
[bdir,bname] = fileparts(Moving_dFile);
if contains(bname,'.')
    [~, bname, ~] = fileparts(bname);
end
tdir = fullfile(bdir,['elxreg_',bname]);
if exist(tdir,'dir')==0
    mkdir(tdir);
end

%% Load Data and Labels
disp([bname,': Fixed and Label'])
regObj.cmiObj(1).loadImg(0,{Fixed_dFile; Fixed_mFile});

disp([bname,': Moving and Label'])
regObj.cmiObj(2).loadImg(0,{Moving_dFile; Moving_mFile});
pause(1);


%% Generate Registration Masks
mask = ismember(regObj.cmiObj(1).img.mat(:,:,:,2),1:60); % create fixed mask
regObj.cmiObj(1).img.mask.merge('replace',mask); % set cmiObj mask

mask = ismember(regObj.cmiObj(2).img.mat(:,:,:,2),1:60);   % Create moving mask for Initial transform
regObj.cmiObj(2).img.mask.merge('replace',mask);


%% Set Parameters
% First make sure the schedule is clear:
regObj.UDschedule;

% Is the fixed image incremental or continuous?
nslc = regObj.cmiObj(1).img.dims(3);
f_sched = {8*ones(1,3), 10*ones(1,3), [5*ones(1,3),2*ones(1,3)]};
if nnz(std(reshape(regObj.cmiObj(1).img.mat(:,:,:,1),[],nslc),[],1)) < nslc
    f_sched{1}(3) = 1;
    f_sched{2}(3) = 1;
    f_sched{3}([3,6]) = 1;
    disp('Fixed image is incremental')
end
nslc = regObj.cmiObj(2).img.dims(3);
m_sched = {8*ones(1,3), 10*ones(1,3), [5*ones(1,3),2*ones(1,3)]};
if nnz(std(reshape(regObj.cmiObj(2).img.mat(:,:,:,1),[],nslc),[],1)) < nslc
    m_sched{1}(3) = 1;
    m_sched{2}(3) = 1;
    m_sched{3}([3,6]) = 1;
    disp('Moving image is incremental')
end

% Add steps to the schedule
regObj.setOdir(tdir);
regObj.addElxStep('Affine','NumberOfResolutions',1,...
                           'FixedImagePyramidSchedule',f_sched{1},...
                           'MovingImagePyramidSchedule',f_sched{1},...
                           'MaximumNumberOfIterations',4000,...
                           'SP_A',50,...
                           'SP_a',2000,...
                           'NumberOfSpatialSamples',2000);
                      
regObj.addElxStep('Warp','NumberOfResolutions',1,...
                         'FixedImagePyramidSchedule',f_sched{2},...
                         'MovingImagePyramidSchedule',m_sched{2},...
                         'MaximumNumberOfIterations',20000,...
                         'SP_A',50,...
                         'SP_a',200000,...
                         'NumberOfSpatialSamples',2000,...
                         'GridSpacingSchedule',10*ones(1,3),...
                         'TransformBendingEnergy',50);
                    
regObj.addElxStep('Warp','NumberOfResolutions',2,...
                         'FixedImagePyramidSchedule',f_sched{3},...
                         'MovingImagePyramidSchedule',m_sched{3},...
                         'MaximumNumberOfIterations',[10000,2000],...
                         'SP_A',[50,50],...
                         'SP_a',[100000,25000],...
                         'NumberOfSpatialSamples',[5000,5000],...
                         'GridSpacingSchedule',[5*ones(1,3),2*ones(1,3)],...
                         'TransformBendingEnergy',25);

%% Run Reg Function then queue
disp([bname,': Initial transform:'])
regObj.VOI2Tform;
pause(1);

disp([bname,': Remove Moving Mask'])
regObj.cmiObj(2).clearMask;

disp([bname,': Enqueue Registration:'])
regObj.runElx(true);
pause(1);



