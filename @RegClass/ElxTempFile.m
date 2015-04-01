function [expfn, insfn, origfn, voifn] = ElxTempFile(self, varargin)
%makes the temporary file directory on one run of Elastix

if self.cmiObj(1).img.check && self.cmiObj(2).img.check
    outdir = self.elObj.wdir;
    ext = '.mhd';
    origfn = fullfile(outdir,['elxtemp-mOrig',ext]);
    insfn = fullfile(outdir,['elxtemp-mImg',ext]);
    expfn = fullfile(outdir,['elxtemp-fImg',ext]);
    voifn = fullfile(outdir,['elxtemp-fMask',ext]);
    if numel(varargin) == 1
        clamp = varargin{1};
    else
        clamp  = 0;
    end
    

    fovR = self.cmiObj(1).img.dims(1:3) .* self.cmiObj(1).img.voxsz;
    fovH = self.cmiObj(2).img.dims(1:3) .* self.cmiObj(2).img.voxsz;    

    vec = self.cmiObj(1).vec;
    
    fnameR = 'elxtemp-fImg';
    fnameH = 'elxtemp-mImg';
    if self.cmiObj(1).img.mask.check
        
        %Set scale if not already there
        scl = -1000*[(min(self.cmiObj(1).img.getMaskVals(self.cmiObj(1).vec))>=0),...
              (min(self.cmiObj(2).img.mat(self.cmiObj(2).img.mat(:,:,:,1)>=self.cmiObj(2).img.valExt(1)))>=0)];
        
        % Un-clamped Moving
        str = cmi_save(0,self.cmiObj(2).img.mat(:,:,:,self.cmiObj(2).vec)+scl(2),...
                        {'Ins'},fovH,origfn);
        disp(['Saved ',str]);
        % Clamped Moving
        timg = self.cmiObj(2).img.mat(:,:,:,self.cmiObj(2).vec)+scl(2); 
        timg(timg>clamp) = 0;
        cmi_save(0,timg,{'Ins'},fovH,insfn);
        disp(['Saved ',str]);
        % Clamped Fixed
        timg = self.cmiObj(1).img.mat(:,:,:,self.cmiObj(1).vec)+scl(1); timg(timg>0) = 0;
        str = cmi_save(0,timg,{'Exp'},fovR,expfn);
        disp(['Saved ',str]);
        % VOI
        str = cmi_save(1,self.cmiObj(1).img.mask.mat(:,:,:,1),{'VOI'},fovR,voifn);
        disp(['Saved ',str]); 
        
    else
        saveMHD(fullfile(self.elObj.wdir, fnameR), self.cmiObj(1).img.mat(:,:,:), self.cmiObj(1).img.labels(vec), fovR);
        saveMHD(fullfile(self.elObj.wdir, fnameH), self.cmiObj(2).img.mat(:,:,:), self.cmiObj(2).img.labels(vec), fovH);
    end
    
end

end

