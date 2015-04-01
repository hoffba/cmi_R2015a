% BH script
function stat = regMouseLungs(cmiObj0,cmiObj1)
%   Input: cmiObj0 = CMIclass object containing Expiration image
%          cmiObj1 = CMIclass object containing Inspiration image

stat = false;
if (nargin==2) && isa(cmiObj0,'CMIclass') && isa(cmiObj1,'CMIclass') ...
        && cmiObj0.img.check && cmiObj1.img.check && cmiObj0.img.mask.check
    
    % Determine filenames
    outdir = cmiObj0.img.dir;
    [~,dstr] = fileparts(outdir);
    if any(strcmpi(dstr(end-2:end),{'exp','ins'}))
        outdir = fileparts(outdir);
    end
    outdir = fullfile(outdir,'elxregRR');
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end
    
    % Temporary files saved as MHD for Elastix compatibility
    ext = '.mhd';
    origfn = fullfile(outdir,['elxtemp-mOrig',ext]);
    insfn = fullfile(outdir,['elxtemp-mImg',ext]);
    expfn = fullfile(outdir,['elxtemp-fImg',ext]);
    voifn = fullfile(outdir,['elxtemp-fMask',ext]);

    % Determine whether images need to be scaled to actual HU values:
    scl = -1000*[(min(cmiObj0.img.getMaskVals(cmiObj0.vec))>=0),...
                 (min(cmiObj1.img.mat(cmiObj1.img.mat(:,:,:,1)>=cmiObj1.img.valExt(1)))>=0)];
    
    % Temporarily save .mhd images
    fov0 = cmiObj0.img.voxsz.*cmiObj0.img.dims(1:3);
    fov1 = cmiObj1.img.voxsz.*cmiObj1.img.dims(1:3);
    disp('Saving temporary image files ...');
    % Un-clamped Ins
    str = cmi_save(0,cmiObj1.img.mat(:,:,:,cmiObj1.vec)+scl(2),...
                    {'Ins'},fov1,origfn);
    disp(['Saved ',str]);
    % Clamped Ins
    timg = cmiObj1.img.mat(:,:,:,cmiObj1.vec)+scl(2); timg(timg>0) = 0;
    str = cmi_save(0,timg,{'Ins'},fov1,insfn);
    disp(['Saved ',str]);
    % Clamped Exp
    timg = cmiObj0.img.mat(:,:,:,cmiObj0.vec)+scl(1); timg(timg>0) = 0;
    str = cmi_save(0,timg,{'Exp'},fov0,expfn);
    disp(['Saved ',str]);
    % VOI
    str = cmi_save(1,cmiObj0.img.mask.mat(:,:,:,1),{'VOI'},fov0,voifn);
    disp(['Saved ',str]);

    % Create Elastix Parameter File:
    p = initElxParStruct; % Initial settings are for mice
    if any(fov0>35) % Human data
        p.Metric.Metric1Weight = 50;
        p.Optimizer.MaximumStepLength = [20,15,15,10];
    end
    pratio = [0.033,0.017,0.008,0.002];
    fips = round(cmiObj0.img.dims(1:3)'*pratio); fips(fips<1)=1;
    mips = round(cmiObj1.img.dims(1:3)'*pratio); mips(mips<1)=1;
    p.Pyramid.FixedImagePyramidSchedule = fips(:);
    p.Pyramid.MovingImagePyramidSchedule = mips(:);
    pstr = saveElxParFile(p,fullfile(outdir,'ElastixParameters.txt'));
    pstraffine = fullfile(fileparts(which('cmi')),'elastix/MouseLungs/MouseLungs.affine.txt');

    % Window Name:
    namestr = ['Elastix Registration: ',cmiObj1.img.name,' --> ',cmiObj0.img.name];
    % Call to Elastix:
    elxstr = ['/opt/elastix/bin/elastix',...
              ' -f ',expfn,...
              ' -m ',insfn,...
              ' -fMask ',voifn,...
              ' -out ',outdir,...
              ' -p ',pstraffine,...
              ' -p ',pstr,'; '];
    % Call to Transformix:
    tfxstr = ['/opt/elastix/bin/transformix -jac all',...
              ' -out ',outdir,...
              ' -in ',origfn,...
              ' -tp ',fullfile(outdir,'TransformParameters.1.txt'),'; '];
    % Cleanup string (removes temporary image files):
    custr = ['find ',outdir,' -name "elxtemp-*" -exec rm -f {} \; '];
    
    % Executes in new xterm window
    %  -geometry sets window size
    %  -T sets window name
    %  -e executes commands (within ' '), must be last xterm input
    %  csh at end to keep window open
    %  & at end to free Matlab from waiting
    stat = system(['xterm -geometry 170x50 -T "',namestr,'"',...
                   ' -e ''',elxstr,tfxstr,custr,';csh''&']);
    disp('check')
end
