function S = imageMeta(fpath)
% Return image metadata for use in MiTAP pipeline

if ischar(fpath)
    if isfolder(fpath)
        % DICOM
    elseif isfile(fpath)
        % Image file
        if endsWith(fpath,'.nii.gz')
            info = niftiinfo(fpath);
        elseif endsWith(fpath,'.mhd')
            info = readMHDinfo(fpath);
        else
            warning('File type ''%s'' not yet programmed.',flip(extractBefore(flip(fpath),'.')));
        end
    else
        warning('Invalid image path: %s',fpath);
    end
else
    error('Invalid input.');
end