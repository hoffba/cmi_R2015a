function stat = regMouseLungElx_batch(muschk,expfn,insfn,voifn,origfn,outdir)

if muschk
    pstr = '/mnt/cmi/projects/CMI/CMI_Matlab_Programs/BenHoff/cmi/elastix/MouseLungs/MouseLungs.bs8-mm.txt';
else
    pstr = '/mnt/cmi/projects/CMI/CMI_Matlab_Programs/BenHoff/cmi/elastix/HumanLungs/HumanLungs.bs8-mm.txt';
end

% Generate call to Elastix:
sysstr = ['/opt/elastix/bin/elastix',...
          ' -f ',expfn,...
          ' -m ',insfn,...
          ' -fMask ',voifn,...
          ' -out ',outdir,...
          ' -p ',pstr];
stat = system(sysstr);

if stat==0
    % Generate call to Transformix for final transformation and |Jac| calculation:
    stat = system(['/opt/elastix/bin/transformix -jac all -out ',outdir,...
                   ' -in ',origfn,...
                   ' -tp ',fullfile(outdir,'TransformParameters.0.txt')]);
    if stat==0
        % Delete temporary files:
        fnames = dir(fullfile(outdir,'*'));
            % **** These are the files to keep:
        fnames(1:2) = []; % '.' and '..'
        fnames(strcmp('TransformParameters.0.txt',{fnames(:).name})) = [];
        fnames(strcmp('elastix.log',{fnames(:).name})) = [];
        fnames(strcmp('transformix.log',{fnames(:).name})) = [];
        fnames(strcmp('result.mhd',{fnames(:).name})) = [];
        fnames(strcmp('result.raw',{fnames(:).name})) = [];
        fnames(strcmp('spatialJacobian.mhd',{fnames(:).name})) = [];
        fnames(strcmp('spatialJacobian.raw',{fnames(:).name})) = [];
            % ****
        fnames = cellfun(@(x)fullfile(outdir,x),{fnames(:).name},'UniformOutput',false);
        delete(fnames{:});
    end
end