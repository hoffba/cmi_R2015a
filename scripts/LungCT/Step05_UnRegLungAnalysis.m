function [expData, insData] = Step05_UnRegLungAnalysis(procdir, fname_ScatNet, regObj)

% calculate densitometry and ScatterNet for whole-lung.
% Note: need to custimize this for lobes.

% generate lung analysis on unreg ins and exp CT scans
% Metrics include: LAA950, HAA810, LAA856, ScatterNet
% idx_Exp = idx_EIV(1);
% idx_Exp_voi = idx_EIV(2);
% idx_Ins = idx_EIV(3);
% idx_Ins_voi = idx_EIV(4);

if ~isfolder(procdir)
    mkdir(procdir);
end

%% ScatterNet for AT on Exp CT scan
SNfname = fullfile(procdir,[fname_ScatNet,'.nii.gz']);
if exist(fullfile(procdir,[fname_ScatNet,'.nii.gz']),'file')
    regObj.cmiObj(1).loadImg(1,SNfname);
%     atMap = niftiread(SNfname);
else
    [atMap, laaMap, percentatMap, percentlaaMap, meanvalatMap, meanvallaaMap] = ScatNet( ...
        regObj.cmiObj(1).img.mat, regObj.cmiObj(1).img.mask.mat, 0);
    ind = regObj.cmiObj(1).imgAppend(atMap,{'ScatNet'});
    regObj.cmiObj(1).img.saveImg(ind,SNfname);
%     niftiwrite(int8(atMap),fullfile(procdir,fname_ScatNet),data(1).voi.info,'Compressed',true);
end

%% Calculate  AT, Emph and PD
filt_flag = length(find(std(single(regObj.cmiObj(1).img.mat(:,:,:,1)),1,[1 2])~=0)') == regObj.cmiObj(1).img.dims(3);

%% First handle Exp data:
fprintf('Analyzing Exp ... ');
% filter data using median filter
temp_data = double(regObj.cmiObj(1).img.mat(:,:,:,1));
if filt_flag
    fprintf('3D Median Filter\n');
    temp_data = medfilt3(temp_data);
else
    fprintf('2D Median Filter for Incremental Scan\n')
    for i_slice = 1: size(temp_data,3)
        temp_data(:,:,i_slice) = medfilt2(temp_data(:,:,i_slice),[3 3]);
    end
end

% Generate mask for whole lung
temp_mask = regObj.cmiObj(1).img.mask.mat;

% Exp CT scan: Vol, Mean, LAA856, SNperc, SNmean
nvox = nnz(temp_mask);
maskvals = temp_data(temp_mask);
expData(1,1) = nvox * prod(regObj.cmiObj(1).img.voxsz)./1e6;
expData(1,2) = mean(maskvals);
expData(1,3) = 100 * nnz(maskvals<-856) / nvox;
expData(1,4) = 100 * nnz(atMap(temp_mask)) / nvox;
expData(1,5) = mean(temp_data(logical(atMap)));


%% Next, handle Ins
fprintf('Analyzing Ins ... ');
% filter data using median filter
temp_data = double(regObj.cmiObj(2).img.mat(:,:,:,1));
if filt_flag
    fprintf('3D Median Filter\n');
    temp_data = medfilt3(temp_data);
else
    fprintf('2D Median Filter for Incremental Scan\n')
    for i_slice = 1: size(temp_data,3)
        temp_data(:,:,i_slice) = medfilt2(temp_data(:,:,i_slice),[3 3]);
    end
end

% Generate masks for whole lung
temp_mask = regObj.cmiObj(2).img.mask.mat;

% Ins CT scan: Vol, Mean, LAA950, LAA810
nvox = nnz(temp_mask);
maskvals = temp_data(temp_mask);
insData(1,1) = nvox * prod(regObj.cmiObj(2).img.voxsz) / 1e6;
insData(1,2) = mean(maskvals);
insData(1,3) = 100 * nnz(maskvals < -950) / nvox;
insData(1,4) = 100 * nnz((maskvals >= -810) & (maskvals < -250)) / nvox;

