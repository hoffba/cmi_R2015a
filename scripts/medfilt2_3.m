function img = medfilt2_3(img)
for i = 1:size(img,3)
    img(:,:,i) = medfilt2(img(:,:,i));
end