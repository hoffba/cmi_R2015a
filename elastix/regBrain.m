% Perform Elastix coregistration with some pre-processing
function stat = regBrain(cmiObj0,cmiObj1,varargin)
%   Input: cmiObj0 = CMIclass object containing fixed image
%          cmiObj1 = CMIclass object containing moving image
%          (Optional) Name/Value pairs:
%               'Affine' = true/false
%               'Debug'  = true/false
%               'Dilate'        = [r1,r2,r3]
%               'Filter'        = [d1,d2,d3], 3D median
%               'Ftype'         = 'median','gauss','average'
%               'Clamp'         = [min,max], truncate image values
%               'T0'            = [s1,s2,s3,t1,t2,t3], initial transform
%               'SP_a'          = [1xR] : R = number of resolutions
%               'SP_A'          = [1xR]
%               'Iterations'    = [1xR]
%               'FixedPyramid'  = [1x(3*R)] (3 values for each)
%               'MovingPyramid' = [1x(3*R)] (3 values for each)

stat = false;
if (nargin>=2) && isa(cmiObj0,'CMIclass') && isa(cmiObj1,'CMIclass') ...
        && cmiObj0.img.check && cmiObj1.img.check
    
    % Parse variable inputs:
    p = inputParser;
    p.CaseSensitive = true;
    addParamValue(p,'Affine',false,@(x)islogical(x)&&(length(x)==1));
    addParamValue(p,'Debug',false,@(x)islogical(x)&&(length(x)==1));
    addParamValue(p,'SP_a',[],@isnumeric);
    addParamValue(p,'SP_A',[],@isnumeric);
    addParamValue(p,'Iterations',[],@isnumeric);
    addParamValue(p,'Dilate',[],@(x)isnumeric(x)&&(ismember(length(x),[1,3])));
    addParamValue(p,'Filter',[],@(x)isnumeric(x)&&(ismember(length(x),[1,3])));
    addParamValue(p,'Ftype','median',@(x)ischar(x)&&any(strcmpi(x,{'median','gauss','average'})));
    addParamValue(p,'Clamp',[],@(x)isnumeric(x)&&(length(x)==2));
    addParamValue(p,'T0',[],@(x)isnumeric(x)&&(length(x)==6));
    addParamValue(p,'FixedPyramid',[],@(x)isnumeric(x)&&(mod(length(x),3)==0));
    addParamValue(p,'MovingPyramid',[],@(x)isnumeric(x)&&(mod(length(x),3)==0));
    parse(p,varargin{:});
    vin = p.Results;
    
    % Determine filenames
    outdir = cmiObj0.img.dir;
    [~,dstr] = fileparts(outdir);
    if any(strcmpi(dstr(end-2:end),{'exp','ins'}))
        % In DICOM directory, so need to look up one level
        outdir = fileparts(outdir);
    end
    outdir = fullfile(outdir,['elxreg_',cmiObj1.img.name]);
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end
    
    % Temporary files saved as MHD for Elastix compatibility
    ext = '.mhd';
    origfn = fullfile(outdir,['elxtemp-mOrig',ext]);
    insfn = fullfile(outdir,['elxtemp-mImg',ext]);
    expfn = fullfile(outdir,['elxtemp-fImg',ext]);
    fmaskfn = '';
    if cmiObj0.img.mask.check
        fmaskfn = fullfile(outdir,['elxtemp-fMask',ext]);
    end
    mmaskfn = '';
    if cmiObj1.img.mask.check
        mmaskfn = fullfile(outdir,['elxtemp-mMask',ext]);
    end

    % Determine whether images need to be scaled to actual HU values:
    scl = -1000*[(min(cmiObj0.img.getMaskVals(cmiObj0.vec))>=0),...
                 (min(cmiObj1.img.mat(cmiObj1.img.mat(:,:,:,1)>=cmiObj1.img.valExt(1)))>=0)];
                 
    % Temporarily save .mhd images
    fov0 = cmiObj0.img.voxsz.*cmiObj0.img.dims(1:3);
    fov1 = cmiObj1.img.voxsz.*cmiObj1.img.dims(1:3);
    disp('Saving temporary image files ...');
% Original moving image
    str = cmi_save(0,cmiObj1.img.mat(:,:,:,cmiObj1.vec)+scl(2),...
                    {'Ins'},fov1,origfn);
    disp(['Saved ',str]);
% Clamped moving image
    timg = cmiObj1.img.mat(:,:,:,cmiObj1.vec);
    if ~isempty(vin.Clamp)
        timg = timg + scl(2);
        timg(timg<vin.Clamp(1)) = vin.Clamp(1);
        timg(timg>vin.Clamp(2)) = vin.Clamp(2);
    end
    str = cmi_save(0,filtImg(timg,vin.Filter,vin.Ftype),{'Ins'},fov1,insfn);
    disp(['Saved ',str]);
% Clamped fixed image
    timg = cmiObj0.img.mat(:,:,:,cmiObj0.vec);
    if ~isempty(vin.Clamp)
        timg = timg+scl(1);
        timg(timg<vin.Clamp(1)) = vin.Clamp(1);
        timg(timg>vin.Clamp(2)) = vin.Clamp(2);
    end
    str = cmi_save(0,filtImg(timg,vin.Filter,vin.Ftype),{'Exp'},fov0,expfn);
    disp(['Saved ',str]);
% Fixed VOI
    if ~isempty(fmaskfn)
        mask = cmiObj0.img.mask.mat(:,:,:,1);
        if ~isempty(vin.Dilate) && all(vin.Dilate>1)
            mask = imdilate(mask,bwellipsoid(vin.Dilate));
        end
        str = cmi_save(1,mask,{'VOI'},fov0,fmaskfn);
        disp(['Saved ',str]);
    end
% Moving VOI
    if ~isempty(mmaskfn)
        mask = cmiObj1.img.mask.mat(:,:,:,1);
        if ~isempty(vin.Dilate) && all(vin.Dilate>1)
            mask = imdilate(mask,bwellipsoid(vin.Dilate));
        end
        str = cmi_save(1,mask,{'VOI'},fov1,mmaskfn);
        disp(['Saved ',str]);
    end

    % Create initial transformation based on VOI extents
    pars = [1 0 0 0 1 0 0 0 1 0 0 0];
    if ~isempty(vin.T0)
        pars([1,5,9:12]) = vin.T0;
    elseif cmiObj0.img.mask.check && cmiObj1.img.mask.check
        % Scale / Translate based on VOI limits
        [v0x,v0y,v0z] = ind2sub(cmiObj0.img.dims(1:3),find(cmiObj0.img.mask.mat));
        [v1x,v1y,v1z] = ind2sub(cmiObj1.img.dims(1:3),find(cmiObj1.img.mask.mat));
        v0x = [min(v0x),max(v0x)] * cmiObj0.img.voxsz(1) - fov0(1)/2;
        v0y = [min(v0y),max(v0y)] * cmiObj0.img.voxsz(2) - fov0(2)/2;
        v0z = [min(v0z),max(v0z)] * cmiObj0.img.voxsz(3) - fov0(3)/2;
        v1x = [min(v1x),max(v1x)] * cmiObj1.img.voxsz(1) - fov1(1)/2;
        v1y = [min(v1y),max(v1y)] * cmiObj1.img.voxsz(2) - fov1(2)/2;
        v1z = [min(v1z),max(v1z)] * cmiObj1.img.voxsz(3) - fov1(3)/2;
        pars([1,5,9]) = [ diff(v1x)/diff(v0x) ,...
                          diff(v1y)/diff(v0y) ,...
                          diff(v1z)/diff(v0z) ];
        pars(10:12) =   [ (sum(v1x)-sum(v0x)*pars(1))/2 ,...
                          (sum(v1y)-sum(v0y)*pars(5))/2 ,...
                          (sum(v1z)-sum(v0z)*pars(9))/2 ];
    else
        % Initialize to 10% x 10% x 40% stretch in ventilating lungs
        pars([1,5,9]) = [1.1,1.1,1.4];
    end
    info = struct('voxsz',cmiObj0.img.voxsz,'dims',cmiObj0.img.dims(1:3));
    inittpfn = fullfile(outdir,'InitialTransform.txt');
    saveAffineTxformPars(pars,info,inittpfn);
    
    % Use hard parameter files:
    elxdir = fullfile(fileparts(which('cmi')),'elastix');
    clinchk = any(fov0>35);
    pstraffine = '';
    if vin.Affine
        if clinchk % Human
            pstraffine = fullfile(elxdir,'HumanLungs','HumanLungs.affine.txt');
    %         pstr = fullfile(elxdir,'HumanLungs','HumanLungs.bs8-mm.txt');
        else % Mouse
            pstraffine = fullfile(elxdir,'MouseLungs','MouseLungs.affine.txt');
    %         pstr = fullfile(elxdir,'MouseLungs','MouseLungs.bs8.txt');
        end
    end
    
    % Create Elastix Parameter File:
%     vargin = {};
%     if ~isempty(vin.Iterations)
%         vargin = {'Iterations',vin.Iterations;
%     end
%     if ~isempty(vin.FixedPyramid)
%         if mod(vin.FixedPyramid,3)==0
%             fR = reshape(repmat(vin.FixedPyramid(:)',3,1),...
%                 1,length(vin.FixedPyramid)*3);
%         else
%             fR = vin.FixedPyramid;
%         end
%         if isempty(nres)
%             nres = length(fR)/3;
%         else
%             nres = min(nres,length(fR));
%         end
%     end
%     if ~isempty(vin.MovingPyramid)
%         if mod(vin.MovingPyramid,3)==0
%             mR = repmat(vin.MovingPyramid(:)',3,1);
%         else
%             mR = vin.MovingPyramid;
%         end
%         if isempty(nres)
%             nres = length(mR)/3;
%         else
%             nres = min(nres,length(mR));
%         end
%     end
%     if ~isempty(vin.SP_A)
%         SP_A = vin.SP_A;
%         if isempty(nres)
%             nres = length(SP_A)/3;
%         else
%             nres = min(nres,length(SP_A));
%         end
%     end
%     if ~isempty(vin.SP_a)
%         SP_a = vin.SP_a;
%         if isempty(nres)
%             nres = length(SP_a)/3;
%         else
%             nres = min(nres,length(SP_a));
%         end
%     end
    p = initElxParStruct(cmiObj0.img.dims(1:3),...
                         cmiObj1.img.dims(1:3),...
                         clinchk,...
                         varargin{:});
    if vin.Debug
        p.Outputs.WriteResultImageAfterEachResolution = 'true';
    end
    
    
    if ~isempty(vin.SP_a)
        nres = min(p.Pyramid.NumberOfResolutions,length(vin.SP_a));
        dnres = p.Pyramid.NumberOfResolutions-nres;
        p.Pyramid.NumberOfResolutions = nres;
        p.Pyramid.FixedImagePyramidSchedule = p.Pyramid.FixedImagePyramidSchedule(3*dnres+1:end);
        p.Pyramid.MovingImagePyramidSchedule = p.Pyramid.MovingImagePyramidSchedule(3*dnres+1:end);
        p.Transform.GridSpacingSchedule = p.Transform.GridSpacingSchedule(3*dnres+1:end);
        p.ImageSampler.NumberOfSpatialSamples = p.ImageSampler.NumberOfSpatialSamples(dnres+1:end);
        p.Optimizer.SP_a = vin.SP_a;
    end
    pstr = saveElxParFile(p,fullfile(outdir,'ElastixParameters.txt'));

    % Window Name:
    namestr = ['Elastix Registration: ',cmiObj1.img.name,' --> ',cmiObj0.img.name];
    % Call to Elastix:
    elxstr = ['/opt/elastix/bin/elastix',...
              ' -f ',expfn,...
              ' -m ',insfn,...
              ' -out ',outdir,...
              ' -t0 ',inittpfn];%,...' -threads 2'];
    if ~isempty(fmaskfn)
        elxstr = [elxstr,' -fMask ',fmaskfn];
    end
    if ~isempty(mmaskfn)
        elxstr = [elxstr,' -mMask ',mmaskfn];
    end
    itp = 0;
    if ~isempty(pstraffine)
        elxstr = [elxstr,' -p ',pstraffine];
        itp = 1;
    end
    elxstr = [elxstr,' -p ',pstr,'; '];
    % Call to Transformix:
    tfxstr = ['/opt/elastix/bin/transformix -jac all',...
              ' -out ',outdir,...
              ' -in ',origfn,...
              ' -tp ',fullfile(outdir,['TransformParameters.',num2str(itp),'.txt']),'; '];
    % Cleanup string (removes temporary image files):
    if vin.Debug
        custr = '';
    else
        custr = ['find ',outdir,' -name "elxtemp-*" -exec rm -f {} \; '];
    end
    
    % Executes in new xterm window
    %  -geometry sets window size
    %  -T sets window name
    %  -e executes commands (within ' '), must be last xterm input
    %  csh at end to keep window open
    %  & at end to free Matlab from waiting
    stat = system(['xterm -geometry 170x50 -T "',namestr,'"',...
                   ' -e ''',elxstr,tfxstr,custr,';csh''&']);
end
end

function img = filtImg(img,n,ftype)
    if ~isempty(n)
        hw = waitbar(0,['Applying ',ftype,' filter, this may take a while ...']);
        switch ftype
            case 'median'
                img = medfilt3(img,n);
            case 'gauss'
            case 'average'
        end
        delete(hw);
    end
end
