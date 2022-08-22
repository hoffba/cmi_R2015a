function result = proc_registration(self)

result = [];

[ct_ref,seg_ref,ct_hom,seg_hom] = self.getData('ct_ref','seg_ref','ct_hom','seg_hom');
if ~(isempty(ct_ref) || isempty(seg_ref) || isempty(ct_hom) || isempty(seg_hom))
    qcheck = false;
    lungreg_BH( ct_ref, self.dat.ct_ref.info, logical(seg_ref),...
                ct_hom, self.dat.ct_hom.info, logical(seg_hom),...
                self.elxdir, self.id, qcheck, self.quickreg);
    % Elastix currently won't save compressed Nifti, so need to compress,
    % rename, and move to procdir
    [~,elxregname] = fileparts(self.getFileName('ct_hom'));
    elxregname = fullfile(self.elxdir,[elxregname,'_R.nii']);
    gzip(elxregname);
    
    % Save registered image to procdir
    movefile([elxregname,'.gz'],self.getFileName('reg_hom'));
end
