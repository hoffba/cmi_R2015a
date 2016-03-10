% RegClass function
function runTfx(self,~,~)
% Run Transformix from existing parameter file

nf = length(self.Tfx.fnames);

if ~exist(self.Tfx.par,'file')
    errordlg('TransformParameter file does not exist.');
elseif ~isdir(self.Tfx.out)
    errordlg('Output directory does not exist.');
elseif (self.Tfx.jac || self.Tfx.jacmat || self.Tfx.def || (nf>0))
    
    self.h.button_TfxSTART.Enable = 'off';
    hw = waitbar(0,'Copying Transform files to new ouput directory');
    [~,pname,ext] = fileparts(self.Tfx.par);
    
    % Save Transforms into output directory
    elxObj = ElxClass;
    elxObj.loadTx(self.Tfx.par);
    % If running on your local iMac, need to adjust filenames for DIPL-run
    % coregistrations:
    if isfield(elxObj.Tx0,'InitialTransformParametersFileName') ...
            && ~strcmp(elxObj.Tx0.InitialTransformParametersFileName,'NoInitialTransform')
        ffname = elxObj.Tx0.InitialTransformParametersFileName;
        if ismac && strncmp(ffname,'/mnt/cmi/',9)
            ffname = ['/Volumes/',ffname(10:end)];
        end
        elxObj.setTx0par('InitialTransformParametersFileName',ffname);
    end
    tpfn = fullfile(self.Tfx.out,[pname,ext]);
    elxObj.saveTx0(tpfn);
    
    tpNN = [self.Tfx.par(1:end-4),'_NN.txt'];
    if any(self.Tfx.nn)
        % Need to save copy of TransformParameter file with NearestNeighbor
        elxObj.setTx0par('ResampleInterpolator','FinalNearestNeighborInterpolator');
        elxObj.saveTx0(tpNN);
    end
    
    waitbar(0.5,hw,'Generating system call to Transformix');
    
    % Start Transformix, run through Matlab since it's fast:
    namestr = '(Transformix) Applying Transform to Selected Images';
    
    ofn = cell(nf,1);
    for i = 1:nf
        [~,tfn,ext] = fileparts(self.Tfx.fnames{i});
        ofn{i} = fullfile(self.Tfx.out,[tfn,'_R',ext]);
    end
    if (nf==0)
        tp = {tpfn};
    else
        tp = cell(nf,1);
        tp(self.Tfx.nn) = {tpNN};
        tp(~self.Tfx.nn) = {tpfn};
    end
    str = elxObj.sysCmd(self.Tfx.out,'title',namestr,'tp',tp,...
        'jac',self.Tfx.jac,'jacmat',self.Tfx.jacmat,'def',self.Tfx.def,...
        'in',self.Tfx.fnames,'outfn',ofn,'wait',true);
    
    % Save command to file for trouble-shooting:
    fid = fopen(fullfile(self.Tfx.out,'transformixCMD.txt'),'w');
    if fid>2
        fprintf(fid,'%s',str);
        fclose(fid);
    end
    % Send command to the system:
    system(str);
    
    delete(hw);
    % CJG added BAH code to correct mhd filename issue.
    fixMHDnames(ofn);
    %--------
    self.h.button_TfxSTART.Enable = 'on';
end
