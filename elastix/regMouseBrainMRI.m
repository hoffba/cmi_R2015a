% BH script
function stat = regMouseBrainMRI(cmiObj0,cmiObj1,waitchk)
%   Input: cmiObj0 = CMIclass object containing Reference image
%          cmiObj1 = CMIclass object containing Homologous image

stat = false;
if (nargin>=2) && isa(cmiObj0,'CMIclass') && isa(cmiObj1,'CMIclass') ...
        && cmiObj0.img.check && cmiObj0.img.mask.check% ...
%         && cmiObj1.img.check && cmiObj1.img.mask.check
    
    waitstr = '''';
    if (nargin<3) || ~waitchk
        waitstr = ';csh''&';
    end
    
    % Determine filenames
    outdir = cmiObj0.img.dir;
%     [~,dstr] = fileparts(outdir);
%     if any(strcmpi(dstr(end-2:end),{'exp','ins'}))
%         outdir = fileparts(outdir);
%     end
    outdir = fullfile(outdir,[cmiObj1.img.name,'_elxreg']);
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end
    
    % ~~~~~~~~ Temporarily save .mhd images
    % Temporary files saved as MHD for Elastix compatibility
    ext = '.mhd';
    orighomfn = fullfile(outdir,['elxtemp-mOrig',ext]);
    homfn = fullfile(outdir,['elxtemp-mImg',ext]);
    reffn = fullfile(outdir,['elxtemp-fImg',ext]);
    fmaskfn = fullfile(outdir,['elxtemp-fMask',ext]);
%     mmaskfn = fullfile(outdir,['elxtemp-mMask',ext]);

    fov0 = cmiObj0.img.voxsz.*cmiObj0.img.dims(1:3);
    fov1 = cmiObj1.img.voxsz.*cmiObj1.img.dims(1:3);
    disp('Saving temporary image files ...');
    
    % Un-processed Homologous image
    str = cmi_save(0,cmiObj1.img.mat(:,:,:,cmiObj1.vec),{'Hom'},fov1,orighomfn);
    disp(['Saved ',str]);
    
    % VOI
    str = cmi_save(1,cmiObj0.img.mask.mat(:,:,:,1),{'VOI'},fov0,fmaskfn);
    disp(['Saved ',str]);
%     str = cmi_save(1,cmiObj1.img.mask.mat(:,:,:,1),{'VOI'},fov1,mmaskfn);
%     disp(['Saved ',str]);
    
    % Processed Homologous image
    cmiObj1.img.imgFilt(1,{cmiObj1.vec*[1,1],[1,3],{'g','u'}});
    % Correct for surface coil signal dropoff
    cmiObj1.img.volCorrectMRI(cmiObj1.vec);
%     timg = cmiObj1.img.mat(:,:,:,cmiObj1.vec); scl = 1000/mean(timg(:));
    str = cmi_save(0,cmiObj1.img.mat(:,:,:,cmiObj1.vec),{'P-Hom'},fov1,homfn);
    disp(['Saved ',str]);
    
    % Processed Reference image
    cmiObj0.img.imgFilt(1,{cmiObj0.vec*[1,1],[1,3],{'g','u'}});
    % Correct for surface coil signal dropoff
    cmiObj0.img.volCorrectMRI(cmiObj0.vec);
%     timg = cmiObj0.img.mat(:,:,:,cmiObj0.vec); scl = 1000/mean(timg(:));
    str = cmi_save(0,cmiObj0.img.mat(:,:,:,cmiObj0.vec),{'P-Ref'},fov0,reffn);
    disp(['Saved ',str]);
    
    % ~~~~~~~~ Generate system call:
    % Create Elastix parameter file:
    p = initAffineParStruct; % Initial settings are for mice
    
    p.Optimizer.Optimizer = 'AdaptiveStochasticGradientDescent';
    p.Optimizer.MaximumStepLength = 0.5;
    p.Optimizer.AutomaticParameterEstimation = 'true';

    p.Pyramid.NumberOfResolutions = 1;
    p.Pyramid.FixedImagePyramidSchedule = [1 1 1];
    p.Pyramid.MovingImagePyramidSchedule = [1 1 1];
    p.Optimizer.SP_a = 5;
    p.Optimizer.SP_alpha = 1;
    p.Optimizer.SP_A = 100;
    
    p.Transform.AutomaticScalesEstimation = 'true';
    p.Transform.Scales = [10 10 10 1 1 1];
    p.Optimizer.MaximumNumberOfIterations = 2000;
    p.ImageSampler.ImageSampler = 'RandomSparseMask';
    p.ImageSampler.NumberOfSpatialSamples = 5000;
    p.Optimizer.NumberOfGradientMeasurements = 10;
%     p.Optimizer.UseAdaptiveStepSizes = 'false';
    pstr = saveElxParFile(p,fullfile(outdir,'ElastixParameters.txt'));
    pstr = fullfile(fileparts(which('cmi')),'elastix/MouseBrainMRI/MouseBrainMRI.euler.txt');

    % Window Name:
    namestr = ['Elastix Registration: ',cmiObj1.img.name,' --> ',cmiObj0.img.name];
    % Call to Elastix:
    elxstr = ['/opt/elastix/bin/elastix',...
              ' -f ',reffn,...
              ' -m ',homfn,...
              ' -fMask ',fmaskfn,...' -mMask ',mmaskfn,...
              ' -out ',outdir,...
              ' -p ',pstr,'; '];
    % Call to Transformix:
    tfxstr = ['/opt/elastix/bin/transformix',...
              ' -out ',outdir,...
              ' -in ',orighomfn,...
              ' -tp ',fullfile(outdir,'TransformParameters.0.txt'),'; '];
    % Cleanup string (removes temporary image files):
    custr = ['find ',outdir,' -name "elxtemp-*" -exec rm -f {} \; '];
    custr = '';
    
    % Executes in new xterm window
    %  -geometry sets window size
    %  -T sets window name
    %  -e executes commands (within ' '), must be last xterm input
    %  csh at end to keep window open
    %  & at end to free Matlab from waiting
    stat = system(['xterm -geometry 170x50 -T "',namestr,'"',...
                   ' -e ''',elxstr,tfxstr,custr,waitstr]);
end
