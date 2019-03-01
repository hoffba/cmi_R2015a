%% Description: Fat-water separation using regularized fieldmap formulation and graph cut solution.
%%
%% Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large
%% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.
%%
%% Some properties:
%%   - Image-space
%%   - 2 species (water-fat)
%%   - Complex-fitting
%%   - Multi-peak fat (pre-calibrated)
%%   - Single-R2*
%%   - Independent water/fat phase
%%   - Requires 3+ echoes at arbitrary echo times (some choices are much better than others! see NSA...)
%%
%%
%% Output: structure outParams
%%   - outParams.species(ii).name: name of the species (taken from algoParams)
%%   - outParams.species(ii).amps: estimated water/fat images, size [nx,ny,ncoils]
%%   - outParams.r2starmap: R2* map (in s^{-1}, size [nx,ny])
%%   - outParams.fieldmap: field map (in Hz, size [nx,ny])
%%
%%
%% Author: Diego Hernando
%% Date created: August 5, 2011
%% Date last modified: November 10, 2011

function [W,F,FF,R2s,dF] = fw_graphcut3d( img, te, f0, df, mask )
% From Diego Hernando's fw_i3cm1i_3pluspoint_hernando_graphcut.m
% Inputs:
%   img = 5D complex image data [nx,ny,nz,nCoils,nTE]
%   te = TE values corresponding to 5th dimension of img
%   f0 = proton frequency
%   df = 3D map of frequency offsets

% Set up parameters:
par = struct('species',{struct('name','water','frequency',0,'relAmps',1),...
                        struct('name','fat','frequency',[3.80, 3.40, 2.60, 1.94, 0.39, -0.60],... % in ppm
                               'relAmps',[0.087 0.693 0.128 0.004 0.039 0.048])},... 
             'size_clique',1,...        % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
             'range_r2star',[0 2500],...% Range of R2* values
             'NUM_R2STARS',101,...      % Number of R2* values for quantization
             'range_fm',[-1000 1000],...% Range of field map values
             'NUM_FMS',101,...          % Number of field map values to discretize
             'NUM_ITERS',40,...         % Number of graph cut iterations
             'SUBSAMPLE',4,...          % Spatial subsampling for field map estimation (for speed)
             'DO_OT',1,...              % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
             'LMAP_POWER',2,...         % Spatially-varying regularization (2 gives ~ uniformn resolution)
             'lambda',0.05,...          % Regularization parameter
             'LMAP_EXTRA',0.05,...      % More smoothing for low-signal regions
             'TRY_PERIODIC_RESIDUAL',0);% Take advantage of periodic residual if uniform TEs (will change range_fm)
[par,d] = validatePar(img,te,df,mask,...
            'size_clique',1,'range_r2star',[0 2500],'range_fm',[-1000 1000],'NUM_FMS',101,'NUM_ITERS',40,...
            'SUBSAMPLE',4,'DO_OT',1,'LMAP_POWER',2,'lambda',0.05,'LMAP_EXTRA',0.05,'TRY_PERIODIC_RESIDUAL',0);
[te,ix] = sort(te);
img = img(:,:,:,:,ix);
        
 % Estimate fieldmap:
 if any(par.SUBSAMPLE>1)
     START = round(par.SUBSAMPLE/2);
     timg = img(START(1):par.SUBSAMPLE(1):d(1),...
                START(2):par.SUBSAMPLE(2):d(2),...
                START(3):par.SUBSAMPLE(3):d(3),:,:);
     tmask = mask(START(1):par.SUBSAMPLE(1):d(1),...
                  START(2):par.SUBSAMPLE(2):d(2),...
                  START(3):par.SUBSAMPLE(3):d(3));
 else
     timg = img;
     tmask = mask;
 end
 
dTE = diff(te);
if par.TRY_PERIODIC_RESIDUAL && sum(abs(dTE-dTE(1)))<1e-6        
    dt = te(2)-te(1);
    period = abs(1/dt);
    NUM_FMS_ORIG = par.NUM_FMS;
    range = diff(par.range_fm);
    par.NUM_FMS = ceil(par.NUM_FMS/range*period);
    par.range_fm = [0 period*(1-1/(par.NUM_FMS))];
    residual = computeResidual( imDataParams, par );
    num_periods = ceil(range/period/2);
    par.NUM_FMS = 2*num_periods*par.NUM_FMS;
    residual = repmat(residual,[2*num_periods 1 1]);
    par.range_fm = [-num_periods*period (num_periods*period-period/NUM_FMS_ORIG)];
else
end






    
    % Compute the residual
    if UNIFORM_TEs % TEST DH* 090801
        % Find out the period in the residual (Assuming uniformly spaced samples)
        dt = imDataParams.TE(2)-imDataParams.TE(1);
        period = abs(1/dt);
        NUM_FMS_ORIG = algoParams.NUM_FMS;
        range = diff(algoParams.range_fm);
        algoParams.NUM_FMS = ceil(algoParams.NUM_FMS/range*period);
        algoParams.range_fm = [0 period*(1-1/(algoParams.NUM_FMS))];
        residual = computeResidual( imDataParams, algoParams );
        num_periods = ceil(range/period/2);
        algoParams.NUM_FMS = 2*num_periods*algoParams.NUM_FMS;
        residual = repmat(residual,[2*num_periods 1 1]);
        algoParams.range_fm = [-num_periods*period (num_periods*period-period/NUM_FMS_ORIG)];
    else
        % If not uniformly spaced TEs, get the residual for the whole range
        residual = computeResidual3d( imDataParams, algoParams );
    end
    

% Setup the estimation, get the lambdamap,...
fms = linspace(algoParams.range_fm(1),algoParams.range_fm(2),algoParams.NUM_FMS);
dfm = fms(2)-fms(1);
lmap = getQuadraticApprox3d( residual, dfm );
lmap = (sqrt(lmap)).^LMAP_POWER;
lmap = lmap + mean(lmap(:))*LMAP_EXTRA;

% Initialize the field map indices
cur_ind = ceil(length(fms)/2)*ones(size(lmap));

% This is the core of the algorithm
fm = graphCutIterations3d(imDataParams,algoParams,residual,lmap,cur_ind );

% If we have subsampled (for speed), let's interpolate the field map
if any(SUBSAMPLE>1)
    fmlowres = fm;
    [SUBY,SUBX,SUBZ] = meshgrid(subY(:),subX(:),subZ(:));
    [ALLY,ALLX,ALLZ] = meshgrid(allY(:),allX(:),allZ(:));
    fm = interpn(SUBX,SUBY,SUBZ,fmlowres,ALLX,ALLY,ALLZ,'spline');
    lmap = interpn(SUBX,SUBY,SUBZ,lmap,ALLX,ALLY,ALLZ,'spline');
    fm(isnan(fm)) = 0;
    imDataParams.images = images0;
    if isfield(imDataParams,'mask')
        imDataParams.mask = mask0;
    end
end
algoParams.lmap = lmap;

% Now take the field map fm and get the rest of the estimates
% w = nan(sx,sy,sz); f = w; r2starmap = w; ff = w;
if algoParams.range_r2star(2)>0
    % DH* 100422 use fine R2* discretization at this point
%     algoParams.NUM_R2STARS = round(algoParams.range_r2star(2)/2)+1;
    outParams.r2starmap = estimateR2starGivenFieldmap3d( imDataParams, algoParams, fm );
else
    outParams.r2starmap = zeros(sx,sy,sz);
end

% if DO_OT % !! NOT adjusted for 3D yet !!
%     imDataParams.fmGC = fm;
% 
%     algoParams.OT_ITERS = 10;
%     algoParams.fieldmap = fm;
%     algoParams.r2starmap = outParams.r2starmap(:,:,:,ka);
%     algoParams.lambdamap = sqrt(lambda*lmap);
% 
%     outParams = fw_i2cm0i_3plusploint_hernando_optimtransfer(imDataParams,algoParams);
%     fm = outParams.fieldmap;
% 
%     % Now re-estimate the R2* map and water/fat images
%     if algoParams.range_r2star(2)>0
%         algoParams.NUM_R2STARS = round(algoParams.range_r2star(2)/2)+1;
%         outParams.r2starmap = estimateR2starGivenFieldmap3d(imDataParams,algoParams,fm);
%     else
%         outParams.r2starmap = zeros(sx,sy,sz);
%     end
% end
amps = decomposeGivenFieldMapAndDampings3d(imDataParams,algoParams,fm,outParams.r2starmap,outParams.r2starmap);
outParams.species(1).amps = amps(:,:,:,1);
outParams.species(2).amps = amps(:,:,:,2);
outParams.fatfract = computeFF(outParams.species(1).amps,outParams.species(2).amps,0);
outParams.fieldmap = fm;
outParams.images = imDataParams.images;
outParams.TE = imDataParams.TE;


function [par,d] = validatePar(img,te,df,mask,varargin)
% 'species',{struct('name','water','frequency',0,'relAmps',1),...
%                         struct('name','fat','frequency',[3.80, 3.40, 2.60, 1.94, 0.39, -0.60],... % in ppm
%                                'relAmps',[0.087 0.693 0.128 0.004 0.039 0.048])},... 
% 'size_clique',1,...        % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
% 'range_r2star',[0 2500],...% Range of R2* values
% 'NUM_R2STARS',101,...      % Number of R2* values for quantization
% 'range_fm',[-1000 1000],...% Range of field map values
% 'NUM_FMS',101,...          % Number of field map values to discretize
% 'NUM_ITERS',40,...         % Number of graph cut iterations
% 'SUBSAMPLE',4,...          % Spatial subsampling for field map estimation (for speed)
% 'DO_OT',1,...              % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
% 'LMAP_POWER',2,...         % Spatially-varying regularization (2 gives ~ uniformn resolution)
% 'lambda',0.05,...          % Regularization parameter
% 'LMAP_EXTRA',0.05,...      % More smoothing for low-signal regions
% 'TRY_PERIODIC_RESIDUAL',0% Take advantage of periodic residual if uniform TEs (will change range_fm)

[d(1),d(2),d(3),d(4),d(5)] = size(img);
fprintf('Image Dimensions: %u %u %u\nNumber of coils: %u\nTE values: %u\n',d(1:3),d(4),d(5));
if numel(te)~=d(5)
    error('Number of TE values does not match image size (5th dimension)');
end

[d2(1),d2(2),d2(3)] = size(df);
if ~isempty(df) && ~all(d(1:3)==d2)
    error('Size of f0 map does not match image');
end

[d2(1),d2(2),d2(3)] = size(mask);
if ~isempty(df) && ~all(d(1:3)==d2)
    error('Size of mask does not match image');
end

p = inputParser;
addParameter(p,'size_clique',   1,              @(x) ismember(x,[1,2]));
addParameter(p,'range_r2star',  [0 2500],       @(x) numel(x)==2 && all(x>=0) && x(1)<=x(2));
addParameter(p,'NUM_R2STARS',   101,            @(x) isnumeric(x) && round(x)>0);
addParameter(p,'range_fm',      [-1000 1000],   @(x) numel(x)==2 && x(1)<=x(2));
addParameter(p,'NUM_FMS',       101,            @(x) isnumeric(x) && round(x)>0);
addParameter(p,'NUM_ITERS',     40,             @(x) isnumeric(x) && round(x)>0);
addParameter(p,'SUBSAMPLE',     4,              @(x) isnumeric(x) && all(round(x)>0) && ismember(numel(x),[1,3]));
addParameter(p,'DO_OT',         1,              @(x) ismember(x,[0,1]));
addParameter(p,'LMAP_POWER',    2,              @(x) ismember(x,[1,2]));
addParameter(p,'lambda',        0.05,           @(x) isnumeric(x) && x>0);
addParameter(p,'LMAP_EXTRA',    0.05,           @(x) isnumeric(x) && x>0);
addParameter(p,'TRY_PERIODIC_RESIDUAL',1,       @(x) ismember(x,[0,1]));
parse(p,varargin{:})
par = p.Results;

par.SUBSAMPLE = round(par.SUBSAMPLE);
if numel(par.SUBSAMPLE)==1
    par.SUBSAMPLE = par.SUBSAMPLE * ones(1,3);
end
par.NUM_R2STARS =  round(par.NUM_R2STARS);
par.NUM_FMS =  round(par.NUM_FMS);
par.NUM_ITERS = round(par.NUM_ITERS);
