function TF = check_CTlung_orientation(img)
% Find orientation of shoulder bones to see if permute is needed

% BAH 2023-02-20                BW = max(img(i).mat(:,:,round(img(i).info.d(3)*4/5):end)>250,[],3);

d = size(img);
se = strel('sphere',3);
BW = max(imclose(img(:,:,round(d(3)*3/4):end)>250,se),[],3);
prop = regionprops(BW,'Orientation','Area');
TF = mod(round(prop([prop.Area]==max([prop.Area])).Orientation/90),2);