function flag = CTlung_ReadImages(ct)

fprintf('\nReading image data from file ... ID = %s\n',ct.id)

%% Check if .exp / .ins files already exist
fnames = fullfile(ct.procdir,strcat(ct.id,'_',ct.tp,'.',{'exp','ins'},ct.ext));
imchk  = cellfun(@(x)exist(x,'file'),fnames);
flag = ~all(imchk);
if flag
    fprintf('   ... from DICOM\n');
    fnames = ct.dcmdir;
else
    fprintf('   ... from NiFTi\n');
end

%% Load images
for i = 1:2
    ct.img(i).loadImg(0,fnames{i});
    
    %% Find orientation of image by bone threshold
    prop = regionprops(max(ct.img(i).mat(:,:,round(img(i).dims(3)/2):end)>800,[],3),...
        'Orientation','Area');
    if mod(round(prop([prop.Area]==max([prop.Area])).Orientation/90),2)
        fprintf('   Permuting image #%u\n',i);
        ct.img(i).permuteMat([2 1 3]);
    end
end

