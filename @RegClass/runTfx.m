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
%         tname = fullfile(self.Tfx.out,'InitialTransform.txt');
%         copyfile(ffname,tname);
        elxObj.setTx0par('InitialTransformParametersFileName',ffname);
    end
    tpfn = fullfile(self.Tfx.out,[pname,ext]);
    elxObj.saveTx0(tpfn);
    
    if any(self.Tfx.nn)
        % Need to save copy of TransformParameter file with NearestNeighbor
        elxObj.setTx0par('ResampleInterpolator','FinalNearestNeighborInterpolator');
        tpNN = [self.Tfx.par(1:end-4),'_NN.txt'];
        elxObj.saveTx0(tpNN);
    end
    
    waitbar(0.5,hw,'Generating system call to Transformix');
    ni = max(1,nf);
    str = cell(ni,1);
    for i = 1:ni
        inC = {};
        tp = tpfn;
        if nf>0
            if self.Tfx.nn(i)
                tp = tpNN;
            end
            [~,oname] = fileparts(self.Tfx.fnames{i});
            inC = {'in',self.Tfx.fnames{i},...
                   'outfn',fullfile(self.Tfx.out,[oname,'_R.mhd'])};
        end
        str(i) = elxObj.tfxCmd(self.Tfx.out,'tp',tp,...
            'jac',self.Tfx.jac,'jacmat',self.Tfx.jacmat,'def',self.Tfx.def,...
            inC{:});
    end
    % Start Transformix, run through Matlab since it's fast:
    namestr = '(Transformix) Applying Transform to Selected Images';
    estr = ['xterm -geometry 170x50 -T "',namestr,'"',...
        ' -e ''',strcat(str{:}),'csh''&'];
    if ismac
        estr = ['/opt/X11/bin/',estr];
    end
    % Save command to file for trouble-shooting:
    fid = fopen(fullfile(self.odir,'transformixCMD.txt'),'w');
    if fid>2
        fprintf(fid,'%s',estr);
        fclose(fid);
    end
    % Send command to the system:
    system(estr);
    
    delete(hw);
    self.h.button_TfxSTART.Enable = 'on';
end
