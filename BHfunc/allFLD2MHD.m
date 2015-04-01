function allFLD2MHD

path = uigetdir(pwd,'Select directory for conversion:');
sfilt = '*.fld';
D = dirtree(path,sfilt);
for i = 1:length(D)
    fnames = dir(fullfile(D{i},sfilt));
    for ifn = 1:length(fnames)
        [~,bname] = fileparts(fnames(ifn).name);
        disp(bname)
        [img,~,fov] = readFLD(fullfile(D{i},fnames(ifn).name));
        if strcmp(bname(end-2:end),'VOI')
            img = logical(img);
        else
            img = img-1000; % this is when back-scaling CT images
        end
        saveMHD(fullfile(D{i},bname),img,{''},fov);
    end
end