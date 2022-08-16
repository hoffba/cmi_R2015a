function [img,swap_flag] = CTlung_readDICOM(dcm_path,tag,sv_path,fn_base,EI_check)
% Output: flag = true: loaded from DICOM ; false: loaded from NIFTI

swap_flag = false;
img = [];

%% Load DICOM image(s)
N = numel(dcm_path);
for i = 1:N
    [img(i).mat,img(i).label,img(i).fov,img(i).orient,img(i).info,~] = cmi_load(1,[],dcm_path{i});

    %% Correct image orientation by bone threshold
    prop = regionprops(max(ct.img(i).mat(:,:,round(img(i).dims(3)/2):end)>800,[],3),'Orientation','Area');
    if mod(round(prop([prop.Area]==max([prop.Area])).Orientation/90),2)
        fprintf('   Permuting image #%u\n',i);
        img(i).mat = permute(img(i).mat,[2 1 3]);
        img(i).fov = img(i).fov([2,1,3]);
    end
end

%% Check EXP/INSP:
if EI_check && numel(img)==2
    swap_flag = nnz(getRespiratoryOrgans(ct.img(1).mat)) > nnz(getRespiratoryOrgans(ct.img(2).mat));
    if swap_flag
        fprintf('Swapping images due to lung volume\n');
        img = img([2,1]);
    end
end

%% Save images as nii
for i = 1:N
    svname = sprintf('%s.%s.nii',fn_base,tag{i});
    fprintf('Saving image to .nii: %s',svname);
    d = size(img(i).img);
    info = init_niftiinfo(sprintf('%s_%s',fn_base,tag{i}),img(i).fov./d,'int16',d);
    niftiwrite(img(1).mat,fullfile(sv_path,svname),info,'Compressed',true);
end


