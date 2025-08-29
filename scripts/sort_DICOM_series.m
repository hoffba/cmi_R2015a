function sort_DICOM_series(dcm_path)

fn = dir(dcm_path);
fn(1:2) = [];
N = numel(fn);

fprintf('Reading DICOM info: %4.1f%%',0);
for i = 1:numel(fn)
    if ~mod(i,round(N/50))
        fprintf('\b\b\b\b\b%4.1f%%',i/N*100);
    end
    if isdicom(fullfile(dcm_path,fn(i).name))
        srcname = fullfile(dcm_path,fn(i).name);
        info = dicominfo(srcname);
        serNum = info.SeriesNumber;
        sdir = fullfile(dcm_path,['series_',num2str(serNum)]);
        if ~isfolder(sdir)
            mkdir(sdir);
        end
        movefile(srcname,sdir);
    end
end
fprintf('\n');

