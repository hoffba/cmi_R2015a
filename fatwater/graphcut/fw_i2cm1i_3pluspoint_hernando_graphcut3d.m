%% Function name: fw_i2cm1i_3pluspoint_hernando
%%
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
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,nz,nTE]
%%   - imDataParams.TE: echo times (in seconds)
%%   - imDataParams.FieldStrength: (in Tesla)
%%
%%   - algoParams.species(ii).name = name of species ii (string)
%%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%%   Example
%%      - algoParams.species(1).name = 'water' % Water
%%      - algoParams.species(1).frequency = [0]
%%      - algoParams.species(1).relAmps = [1]
%%      - algoParams.species(2).name = 'fat' % Fat
%%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
%%
%%   - algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
%%   - algoParams.range_r2star = [0 0]; % Range of R2* values
%%   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
%%   - algoParams.range_fm = [-400 400]; % Range of field map values
%%   - algoParams.NUM_FMS = 301; % Number of field map values to discretize
%%   - algoParams.NUM_ITERS = 40; % Number of graph cut iterations
%%   - algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
%%   - algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
%%   - algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
%%   - algoParams.lambda = 0.05; % Regularization parameter
%%   - algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
%%   - algoParams.TRY_PERIODIC_RESIDUAL = 0; % Take advantage of periodic residual if uniform TEs (will change range_fm)
%%   - algoParams.residual: in case we pre-computed the fit residual (mostly for testing)
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

function outParams = fw_i2cm1i_3pluspoint_hernando_graphcut3d( imDataParams, algoParams )

DEBUG = 0;

% Check validity of params, and set default algorithm parameters if not provided
[validParams,algoParams] = checkParamsAndSetDefaults_graphcut( imDataParams,algoParams );
if validParams==0
    disp('Exiting -- data not processed');
    outParams = [];
    return;
end

% Get data dimensions
[sx,sy,sz,N] = size(imDataParams.images);

% If precession is clockwise (positive fat frequency) simply conjugate data
if imDataParams.PrecessionIsClockwise <= 0
    imDataParams.images = conj(imDataParams.images);
    imDataParams.PrecessionIsClockwise = 1;
end

% Check spatial subsampling option (speedup ~ quadratic SUBSAMPLE parameter)
SUBSAMPLE = ones(1,3);
ns = min(numel(algoParams.SUBSAMPLE),3);
if ns==1
    SUBSAMPLE = algoParams.SUBSAMPLE*SUBSAMPLE;
else
    SUBSAMPLE(1:ns) = algoParams.SUBSAMPLE(1:ns);
end
if any(SUBSAMPLE > 1)
    images0 = imDataParams.images;
    START = round(SUBSAMPLE/2);
    allX = 1:sx;
    allY = 1:sy;
    allZ = 1:sz;
    subX = START:SUBSAMPLE:sx;
    subY = START:SUBSAMPLE:sy;
    subZ = START:SUBSAMPLE:sz;
    imDataParams.images = images0(subX,subY,subZ,:,:);
    if isfield(imDataParams,'mask')
        mask0 = imDataParams.mask;
        imDataParams.mask = imDataParams.mask(subX,subY,subZ,:,:);
    end
end

% Regularization parameter
% lambda = algoParams.lambda;

% Spatially-varying regularization.  The LMAP_POWER applies to the
% sqrt of the curvature of the residual, and LMAP_POWER=2 yields
% approximately uniform resolution.
LMAP_POWER = algoParams.LMAP_POWER;

% LMAP_EXTRA: Extra flexibility for including prior knowledge into
% regularization. For instance, it can be used to add more smoothing
% to noise regions (by adding, eg a constant LMAP_EXTRA), or even to
% add spatially-varying smoothing as a function of distance to
% isocenter...
LMAP_EXTRA = algoParams.LMAP_EXTRA;

% Finish off with some optimization transfer -- to remove discretization
% DO_OT = algoParams.DO_OT;

% Let's get the residual. If it's not already in the params, compute it
if isfield(algoParams,'residual')
    % Grab the residual from the params structure
    residual = algoParams.residual;
else
    % Check for uniform TE spacings
    dTE = diff(imDataParams.TE);
    
    if isfield(algoParams,'TRY_PERIODIC_RESIDUAL')
        TRY_PERIODIC_RESIDUAL = algoParams.TRY_PERIODIC_RESIDUAL;
    else
        TRY_PERIODIC_RESIDUAL=1;
    end
    
    if TRY_PERIODIC_RESIDUAL==1 && sum(abs(dTE - dTE(1)))<1e-6 % If we have uniform TE spacing
        UNIFORM_TEs = 1;
    else
        UNIFORM_TEs = 0;
    end
    
    
    if DEBUG == 1
        fprintf('UNIFORM_TEs = %u',UNIFORM_TEs);
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



