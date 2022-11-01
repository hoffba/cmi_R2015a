function pipeline_reg(img,elxdir,ID,opts_in) 
% Performs lung CT image registration Ins to Exp
%
% Inputs:
%   img : structure containing
%       .flag    = flag for available data
%       .mat     = image matrix
%       .info    = image information
%       .label   = segmentation map
%   elxdir : location for saving elastix files
%   ID : string identifier for image data
%   opts : structure with registration options
%       .fn_log     = log filename
%       .quickreg   = TF for faster third reg step
%       .jac        = TF for saving |Jacobian| map
%       .jacmat     = TF for saving full Jacobian matrix

opts = struct('fn_log','',...
              'quickreg',false,...
              'jac',true,...
              'jacmat',false,...
              'def',false);
if nargin==4
    fn = fieldnames(opts);
    for i = 1:numel(fn)
        if isfield(opts_in,fn{i})
            opts.(fn{i}) = opts_in.(fn{i});
        end
    end
end

%% Set up Reg object:

regObj = RegClass(false);
regObj.cmiObj(1).setImg(img(1).mat,img(1).info.label,img(1).info.fov,img(1).info.orient,img(1).info.name);
regObj.cmiObj(1).img.mask.merge('replace',logical(img(1).label));
regObj.cmiObj(2).setImg(img(2).mat,img(2).info.label,img(2).info.fov,img(2).info.orient,img(2).info.name);
regObj.cmiObj(2).img.mask.merge('replace',logical(img(2).label));

%% Is the fixed image incremental or continuous?
nslc = regObj.cmiObj(1).img.dims(3);
gapchk = nnz(std(reshape(regObj.cmiObj(1).img.mat(:,:,:,1),[],nslc),[],1)) < nslc;
if gapchk
    writeLog(opts.fn_log,'Fixed image is incremental\n');
end
nslc = regObj.cmiObj(2).img.dims(3);
gapchk(2) = nnz(std(reshape(regObj.cmiObj(2).img.mat(:,:,:,1),[],nslc),[],1)) < nslc;
if gapchk(2)
    writeLog(opts.fn_log,'Moving image is incremental\n');
end

%% Set up preprocessing options
regObj.UDpreproc('filtn',[3,3,0;3,3,0],'dilaten',[15,15,15;0,0,0],'clamp',[-1000,0;-1000,0])

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
if opts.quickreg
    nres = 1;
end

if ~isfolder(elxdir)
    mkdir(elxdir)
end
regObj.setOdir(elxdir);

% Clear elastix directory for a fresh start
delete(fullfile(elxdir,'*.*'));

%% Add steps to the schedule
% regObj.addElxStep('Affine','NumberOfResolutions',1,...
%                            'FixedImagePyramidSchedule',f_sched{1},...
%                            'MovingImagePyramidSchedule',f_sched{1},...
%                            'MaximumNumberOfIterations',4000,...
%                            'SP_A',50,...
%                            'SP_a',2000,...
%                            'NumberOfSpatialSamples',2000,...
%                            'DefaultPixelValue',-2000,...
%                            'WriteResultImage','false');
                      
regObj.addElxStep('Warp','NumberOfResolutions',1,...
                         'FixedImagePyramidSchedule',f_sched{2},...
                         'MovingImagePyramidSchedule',m_sched{2},...
                         'MaximumNumberOfIterations',20000,...
                         'SP_A',50,...
                         'SP_a',200000,...
                         'NumberOfSpatialSamples',2000,...
                         'GridSpacingSchedule',10*ones(1,3),...
                         'TransformBendingEnergy',50,...
                         'DefaultPixelValue',-2000,...
                         'WriteResultImage','false');
                    
regObj.addElxStep('Warp','NumberOfResolutions',nres,...
                         'FixedImagePyramidSchedule',f_sched{3}(1:nres),...
                         'MovingImagePyramidSchedule',m_sched{3}(1:nres),...
                         'MaximumNumberOfIterations',maxiter(1:nres),...
                         'SP_A',SP_A(1:nres),...
                         'SP_a',SP_a(1:nres),...
                         'NumberOfSpatialSamples',samp(1:nres),...
                         'GridSpacingSchedule',gridsp(1:nres*3),...
                         'TransformBendingEnergy',25,...
                         'DefaultPixelValue',-2000,...
                         'WriteResultImage','false');

%% Set up other options
opt_str = {'jac','jacmat','def','Tvoi','Tsurf'};
for i = 1:numel(opt_str)
    if isfield(opts,opt_str{i})
        regObj.UDschedule(opt_str{i},opts.(opt_str{i}));
    end
end
                     
%% Run Reg Function then queue
writeLog(opts.fn_log,[ID,': Initial transform...\n']);
regObj.VOI2Tform;
pause(1);

writeLog(opts.fn_log,[ID,': Remove Moving Mask...\n']);
regObj.cmiObj(2).clearMask;

writeLog(opts.fn_log,[ID,': Start Registration...\n']);
regObj.setWait(true);
regObj.setLogFile(opts.fn_log);
[stat,cmdstr] = regObj.runElx(false);
writeLog(opts.fn_log,'Elastix command:\n%s\n',cmdstr);
pause(1);

delete(regObj);



