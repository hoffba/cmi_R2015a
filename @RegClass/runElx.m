% RegClass function
function stat = runElx(self,hObject,~)

if isa(hObject,'matlab.ui.control.UIControl')
    switch hObject.Tag
        case 'button_Start'
            qchk = false;
        case 'button_Queue'
            qchk = true;
    end
elseif isnumeric(hObject)
    qchk = logical(hObject);
elseif islogical(hObject)
    qchk = hObject;
else
    qchk = false;
end

stat = self.cmiObj(1).img.check && self.cmiObj(2).img.check;
hw = waitbar(0,'Setting up Elastix inputs ... Initial Transform');

elxC = {};

% ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
% ThinPlateSpline does not work as initial transform!!
% ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
% % If warping is used for initial transform (ThinPlateSpline), 
% %       points need to be saved in .txt for additional Elastix input
% if ~isempty(self.elxObj.Tx0) && strcmp(self.elxObj.Tx0.Transform,'SplineKernelTransform')
%     fname = self.savePoints(1,fullfile(self.odir,'ipp.txt'));
%     elxC = [elxC,{'ipp'},{fname}];
% end

% Check output directory:
if stat && isempty(self.odir)
    self.setOdir;
end
stat = ~isempty(self.odir);

% Adjust Elastix parameters based on input images:
labl = {'Fixed','Moving'};
npar = length(self.elxObj.Schedule);
for i = 1:2
    if self.cmiObj(i).img.dims(3)<5
        str = 'Recursive';
    else
        str = 'Smoothing';
    end
    self.elxObj.setPar(1:npar,[labl{i},'ImagePyramid'],[labl{i},str,'ImagePyramid']);
end

% Pre-processing & save temporary .mhd files
if stat
    
    % Get FOV for both images:
    ffov = self.cmiObj(1).img.dims(1:3).*self.cmiObj(1).img.voxsz;
    mfov = self.cmiObj(2).img.dims(1:3).*self.cmiObj(2).img.voxsz;
       
    % Set filenames:
    outfn = fullfile(self.odir,[self.cmiObj(2).img.name,'_R.mhd']);
    origfn = fullfile(self.odir,'elxtemp-origm.mhd');
    fname = fullfile(self.odir,'elxtemp-f.mhd');
    mname = fullfile(self.odir,'elxtemp-m.mhd');
    fmskname = fullfile(self.odir,'elxtemp-fMask.mhd');
    mmskname = fullfile(self.odir,'elxtemp-mMask.mhd');
    
    % Filter settings:
    if any(self.filtN>0)
        filtN = 2*self.filtN + 1;
        switch self.ftype
            case 'Gaussian'
                filtN = self.filtN/3;
                func = @(x)imfilter(x,gaussND(filtN));
            case 'Median'
                func = @(x)medfilt2(x,filtN(1:2));
            case 'Wiener'
                func = @(x)wiener2(x,filtN(1:2));
        end
    else
        filtN = [];
    end
    
    % Save original moving image:
    waitbar(0.1,hw,'Saving Original Moving Image ...');
    timg = self.cmiObj(2).img.mat(:,:,:,self.cmiObj(2).vec);
    stat = saveMHD(origfn,timg,[],mfov);
    
    % Moving Image:
    if stat
        if ~isempty(filtN) && any(filtN)
            % Apply image filter
            waitbar(0.25,hw,'Filtering Moving Image ...');
            if strcmp(self.ftype,'Gaussian')
                timg = feval(func,timg);
            else
                for islc = 1:size(timg,3)
                    timg(:,:,islc) = feval(func,timg(:,:,islc));
                    waitbar(islc/size(timg,3),hw);
                end
            end
        end
        % Apply image crops:
        timg(timg>self.clamp(2,2)) = self.clamp(2,2);
        timg(timg<self.clamp(2,1)) = self.clamp(2,1);
        % Apply histogram equalization:
        if self.histeq
            if self.cmiObj(2).img.mask.check
                tmask = self.cmiObj(2).img.mask.mat;
            else
                tmask = [];
            end
            if ~isempty(tmask)
                timg = timg.*tmask;
            end
            mmin = min(timg(:));
            mmax = max(timg(:));
            timg = (timg-mmin)/(mmax-mmin);
            for i = 1:size(timg,3)
                timg(:,:,i) = adapthisteq(timg(:,:,i))*10^4;
                if ~isempty(tmask)
                    timg(:,:,i) = timg(:,:,i).*tmask(:,:,i);
                end
            end
        end
        % Save moving image:
        waitbar(0.35,hw,'Saving Processed Moving Image ...');
        stat = saveMHD(mname,timg,[],mfov);
    end
        
    % Fixed Image:
    if stat
        timg = self.cmiObj(1).img.mat(:,:,:,self.cmiObj(1).vec);
        if ~isempty(filtN) && any(filtN)
            waitbar(0.5,hw,'Filtering Fixed Image ...');
            if strcmp(self.ftype,'Gaussian')
                timg = feval(func,timg);
            else
                for islc = 1:size(timg,3)
                    timg(:,:,islc) = feval(func,timg(:,:,islc));
                    waitbar(islc/size(timg,3),hw);
                end
            end
        end
        % Apply image crops:
        timg(timg>self.clamp(1,2)) = self.clamp(1,2);
        timg(timg<self.clamp(1,1)) = self.clamp(1,1);
        % Apply histogram equalization:
        if self.histeq
            if self.cmiObj(1).img.mask.check
                tmask = self.cmiObj(1).img.mask.mat;
            else
                tmask = [];
            end
            if ~isempty(tmask)
                timg = timg.*tmask;
            end
            mmin = min(timg(:));
            mmax = max(timg(:));
            timg = (timg-mmin)/(mmax-mmin);
            for i = 1:size(timg,3)
                timg(:,:,i) = adapthisteq(timg(:,:,i))*10^4;
                if ~isempty(tmask)
                    timg(:,:,i) = timg(:,:,i).*tmask(:,:,i);
                end
            end
        end
        % Save fixed image:
        waitbar(0.6,hw,'Saving Processed Fixed Image ...');
        stat = saveMHD(fname,timg,[],ffov);
    end
    
    % Dilate and save moving mask:
    if stat
        if self.cmiObj(2).img.mask.check
            waitbar(0.75,hw,['Dilating and Saving Moving VOI ...',mmskname]);
            elxC = [elxC,'mMask',mmskname];
            stat = saveMHD(mmskname,...
                           voiDilate(self.cmiObj(2).img.mask.mat,self.dilateN),...
                           [],mfov);
        end
    end
    
    % Dilate and save fixed mask:
    if stat
        if self.cmiObj(1).img.mask.check
            waitbar(0.85,hw,['Dilating and Saving Moving VOI ...',fmskname]);
            elxC = [elxC,'fMask',fmskname];
            stat = saveMHD(fmskname,...
                           voiDilate(self.cmiObj(1).img.mask.mat,self.dilateN),...
                           [],ffov);
        end
    end
    
end

if stat
    
    % Name of xterm window:
    namestr = ['Elastix Registration: ',self.cmiObj(2).img.name,...
                                ' --> ',self.cmiObj(1).img.name];
    % Final Transformix setup:
    tpfname = fullfile(self.odir,['TransformParameters.',...
        num2str(length(self.elxObj.Schedule)-1),'.txt']);
                    
    % Waiting for completion / running independently:
    waitchk = qchk || self.waitchk;
    
    % Use all cores (max number of threads):
    ncores = feature('numCores');
    
    % Generate system call to elastix/transformix:
    % ~*~ Assumes all relevant files are in directory "odir"
    
    cmdstr = self.elxObj.sysCmd(self.odir,'title',namestr,...
        'f',fname,...
        'm',mname,...
        elxC{:},...
        'in',origfn,...
        'outfn',outfn,...
        'tp',tpfname,...
        'jac',self.jac,'jacmat',self.jacmat,'def',self.def,...
        'cleanup',true,'wait',waitchk,'threads',ncores);
    
    % Save command to file for trouble-shooting:
    fid = fopen(fullfile(self.odir,'elastixCMD.txt'),'w');
    if fid>2
        fprintf(fid,'%s',cmdstr);
        fclose(fid);
    end
    
    % Determine whether to "Start" or "Add to Queue"
    if qchk
        % Add command to Qfile:
        fid = fopen(self.qfile,'at'); % 'at' = append text
        fprintf(fid,'%s\n',cmdstr);
        fclose(fid);
        
        % Remove finished/failed/deleted jobs
        if ~isempty(self.job)
            self.job(~ismember({self.job(:).State},{'running','queued'})) = [];
        end
        
        % Determine if new batch job needs to be started:
        nj = length(self.job);
        if (nj<self.qnum) && exist(self.qfile,'file')
            % Start new batch job:
            tjob = runBatchFromQ(self.qfile);
            if isempty(self.job)
                self.job = tjob;
            else
                self.job(nj+1) = tjob;
            end
        end
        disp(['Added to queue: ',namestr]);
    else
        stat = ~system(cmdstr);
    end
end
delete(hw);
pause(0.01);
end

function mask = voiDilate(mask,r)
    if any(r)
        ni = max(ceil(r/5));
        for j = 1:ni
            rt = min(r,5);
            r = r - rt;
            se = bwellipsoid(rt);
            mask = imdilate(mask,se);
        end
    end
end
