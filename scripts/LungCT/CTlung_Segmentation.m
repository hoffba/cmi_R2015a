function err = CTlung_Segmentation(data,method)

err = [];
try

% Check for GPU (required if method==3)
if method==3 && ~gpuDeviceCount
    warning('DL_lung_segmentation was selected, but no GPU is available. Reverting to getRespiratoryOrgans method.');
    method = 2;
end

for i = 1:numel(imgObj)
    switch method
        case 1
            fprintf('- Generating VOI from Step02_segLungHuman_cjg ...\n');
            tmask = segLungHuman_cjg_bh(1,data.img(i).mat,dname);
        case 2
            fprintf('- Generating VOI from getRespiratoryOrgans ...\n');
            tmask = getRespiratoryOrgans(data.img(i).mat);
        case 3
            fprintf('- Generating VOI from DL_lung_segmetation ...\n');
            tmask = DL_lung_segmetation(data.img(i).mat);
        otherwise
    end
    if ~isempty(tmask)
        imgObj(i).mask.merge('replace',tmask);
    end
end

catch err
end