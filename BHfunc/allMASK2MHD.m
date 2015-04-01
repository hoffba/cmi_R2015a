function allMASK2MHD

path = uigetdir(pwd,'Select directory for conversion:');
sfilt = '*.mask';
D = dirtree(path,sfilt);
for i = 1:length(D)
    fnames = dir(fullfile(D{i},sfilt));
    for ifn = 1:length(fnames)
        [~,bname] = fileparts(fnames(ifn).name);
        disp(bname)
        imgfn = fullfile(D{i},[bname,'.mhd']);
        if strcmp(imgfn(end-7:end-4),'_VOI')
            imgfn(end-7:end-4) = [];
        end
        if exist(imgfn,'file')
            [img,~,fov] = readMHD(imgfn);
            d = size(img);
            mask = readMASK(fullfile(D{i},fnames(ifn).name),d);
            saveMHD(fullfile(D{i},bname),mask,{''},fov);
        else
            disp('-- image file not found!');
        end
    end
end