function result = calcFatMGE(imData)
% Calculate fat fraction map from MGE pulse sequence
% Uses Matlab code from Diego Hernando, a toolbox developed as an initiative of 
%   the ISMRM Fat-Water MRI Workshop

% * Multi-point fat-water separation with R2* using a graphcut field map estimation 
%   (D. Hernando) hernando/graphcut/fw_i2cm1i_3pluspoint_hernando_graphcut.m

if nargin==0
    % User finds file to load:
    %   .mat for loading imData structure
    %   .MRD or fid for loading raw data and reconstructing
    [fname,fdir] = uigetfile('*','Select image data:');
    if ~ischar(fdir)
        result = [];
        return;
    elseif strcmp(fname(end-3:end),'.mat')
        imData = load(fullfile(fdir,fname),'-str');
    elseif any(strcmp(fname(end-2:end),{'fid','MRD'}))
        imData = recon4fw(fullfile(fdir,fname));
    else
        error('Invalid file selected.');
    end
end

if ~all(isfield(imData,{'images','TE','FieldStrength'}))
    error('Invalid image data structure.')
end

% Define algorithm parameters:
spec = struct('name',... name of species ii (string)
                {'water','fat'},... 
              'frequency',... frequency shift in ppm of each peak within species ii
                {0,[3.80, 3.40, 2.60, 1.94, 0.39, -0.60]},... 
              'relAmps',... relative amplitude (sum normalized to 1) of each peak within species ii
                {1,[0.087 0.693 0.128 0.004 0.039 0.048]});
algPar = struct('species',      {spec},... 
                'size_clique',  {1},... Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
                'range_r2star', {[0,0]},... Range of R2* values
                'NUM_R2STARS',  {1},... Numbre of R2* values for quantization
                'range_fm',     {[-400,400]},... Range of field map values
                'NUM_FMS',      {301},... Number of field map values to discretize
                'NUM_ITERS',    {40},... Number of graph cut iterations
                'SUBSAMPLE',    {2},... Spatial subsampling for field map estimation (for speed)
                'DO_OT',        {1},... 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
                'LMAP_POWER',   {2},... Spatially-varying regularization (2 gives ~ uniformn resolution)
                'lambda',       {0.05},... Regularization parameter
                'LMAP_EXTRA',   {0.05},... More smoothing for low-signal regions
                'TRY_PERIODIC_RESIDUAL',{0});% Take advantage of periodic residual if uniform TEs (will change range_fm)  

% Perform Fat-Water analysis:
result = fw_i3cm1i_3pluspoint_hernando_graphcut(imData,algPar);
% Currently set up for single slice, so need a loop:
% tdata = imData;
% for i = 1:size(imData.images,3)
%     tdata.images = imData.images(:,:,i,:,:);
%     result(i) = fw_i2cm1i_3pluspoint_hernando_graphcut(tdata,algPar);
% end


