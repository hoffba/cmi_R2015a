% CMI script
function stat = regLungsElx(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings

stat = false;
rlabels = {'CT-Ins','CT-Exp'};
if cmiObj.img.check && (cmiObj.img.dims(4)==2) && cmiObj.img.mask.check ...
        &&  all(ismember(cmiObj.img.labels,rlabels))
    
    scl = [0,0];
    scl(min(cmiObj.img.getMaskVals,[],1)>=0) = -1000;
    insind = find(strcmp('CT-Ins',cmiObj.img.labels),1);
    expind = find(strcmp('CT-Exp',cmiObj.img.labels),1);
    
    if ~(isempty(insind) || isempty(expind))
        % Determine filenames
        outdir = cmiObj.img.dir;
        [~,dstr] = fileparts(outdir);
        if any(strcmpi(dstr(end-2:end),{'Exp','Ins'}))
            outdir = fileparts(outdir);
        end
        outdir = fullfile(outdir,'elxreg');
        if ~exist(outdir,'dir')
            mkdir(outdir);
        end
        % Temporary files saved as MHD for Elastix compatibility
        origfn = fullfile(outdir,'ct-insOrig.mhd');
        insfn = fullfile(outdir,'ct-ins.mhd');
        expfn = fullfile(outdir,'ct-exp.mhd');
        voifn = fullfile(outdir,'ct-voi.mhd');

        % Temporarily save .mhd images
        fov = cmiObj.img.voxsz.*cmiObj.img.dims(1:3);
        disp('Saving temporary image files ...');
        % Un-clamped Ins
        str = cmi_save(0,cmiObj.img.mat(:,:,:,insind)+scl(insind),{'Ins'},fov,origfn);
        disp(['Saved ',str]);
        % Clamped Ins
        timg = cmiObj.img.mat(:,:,:,insind)+scl(insind); timg(timg>0) = 0;
        str = cmi_save(0,timg,{'Ins'},fov,insfn);
        disp(['Saved ',str]);
        % Clamped Exp
        timg = cmiObj.img.mat(:,:,:,expind)+scl(expind); timg(timg>0) = 0;
        str = cmi_save(0,timg,{'Exp'},fov,expfn);
        disp(['Saved ',str]);
        % VOI
        str = cmi_save(1,cmiObj.img.mask.mat(:,:,:,1),{'VOI'},fov,voifn);
        disp(['Saved ',str]);

        muschk = true;
        if any(cmiObj.img.voxsz>0.1)
            muschk = false; % Human
        end
        j = batch(@regLungsElx_batch,1,{muschk,expfn,insfn,voifn,origfn,outdir});
        disp(['Job (',num2str(j.ID),') started: ',cmiObj.img.dir])
    end
    
%     % Generate call to Elastix:
%     sysstr = ['/opt/elastix/bin/elastix',...
%               ' -f ',expfn,...
%               ' -m ',insfn,...
%               ' -fMask ',voifn,...
%               ' -out ',outdir,...
%               ' -p /mnt/cmi/projects/CMI/CMI_Matlab_Programs/BenHoff/cmi/elastix/MouseLungs/MouseLungs.bs8-mm.txt'];
%     stat = system(sysstr);
%     
%     if stat==0
%         % Generate call to Transformix for final transformation and |Jac| calculation:
%         stat = system(['/opt/elastix/bin/transformix -jac all -out ',fdir,...
%                        ' -in ',origfn,...
%                        ' -tp ',fullfile(outdir,'TransformParameters.0.txt')]);
%         if stat==0
%             stat = rmdir(outdir,'s');
%         end
%     end
    
end
