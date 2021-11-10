function img = CTlung_checkEI(img)

if nnz(regObj.cmiObj(1).img.mask.mat) > nnz(regObj.cmiObj(2).img.mask.mat)
    regObj.swapCMIdata;
end
%% Save nii.gz files using ID and Tag
fprintf('Saving Exp/Ins images to: %s\n',procdir);
for i = 1:2
    regObj.cmiObj(i).img.saveImg(1,fnames{i,1},1);
    regObj.cmiObj(i).img.saveMask(fnames{i,2});
end