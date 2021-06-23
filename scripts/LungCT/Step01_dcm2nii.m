function data = Step01_dcm2nii(procdir, fnames, tags)

% Use dicom catalog csv file to extract information.
% home_pwd = home directory
% ID_list = current ID being processed
% list = files associated with ID (dicom Exp and dicom Ins)

if numel(fnames) ~= numel(tags)
    error('Inputs fnames and tags must have same number of elements')
end
outFormat = '.nii.gz';
data = struct('tag',tags,'img',{[],[]},'voi',{[],[]});
nf = length(fnames);
for i=1:nf
    if ~isfolder(procdir)
        mkdir(procdir)
    end
    
    dicm2nii(fnames{i}, procdir, outFormat)
    list_file = dir(fullfile(procdir,'*.gz'));
    
    if ~isempty(list_file)
        fprintf('Reading image from file (#%u): %s\n',i,list_file.name);
        str_file = fullfile(list_file.folder,list_file.name);

        data(i).img.info = niftiinfo(str_file);
        data(i).img.mat = niftiread(data(i).img.info);
        
        data(i).voi.info = data(i).img.info;
        data(i).voi.info.Datatype = 'int8';
        data(1).voi.info.BitsPerPixel = 8;
        
        delete(str_file)
    end
end