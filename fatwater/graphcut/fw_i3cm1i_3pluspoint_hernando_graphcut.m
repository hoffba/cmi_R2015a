% Function name: fw_i2cm1i_3pluspoint_hernando
%
% Description: Fat-water separation using regularized fieldmap formulation and graph cut solution.
%
% Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large
% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.
%
% Some properties:
%   - Image-space
%   - 2 species (water-fat)
%   - Complex-fitting
%   - Multi-peak fat (pre-calibrated)
%   - Single-R2*
%   - Independent water/fat phase
%   - Requires 3+ echoes at arbitrary echo times (some choices are much better than others! see NSA...)
%
% Input: structures imDataParams and algoParams
%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%   - imDataParams.TE: echo times (in seconds)
%   - imDataParams.FieldStrength: (in Tesla)
%
%   - algoParams.species(ii).name = name of species ii (string)
%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%   Example
%      - algoParams.species(1).name = 'water' % Water
%      - algoParams.species(1).frequency = [0]
%      - algoParams.species(1).relAmps = [1]
%      - algoParams.species(2).name = 'fat' % Fat
%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
%
%   - algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
%   - algoParams.range_r2star = [0 0]; % Range of R2* values
%   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
%   - algoParams.range_fm = [-400 400]; % Range of field map values
%   - algoParams.NUM_FMS = 301; % Number of field map values to discretize
%   - algoParams.NUM_ITERS = 40; % Number of graph cut iterations
%   - algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
%   - algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
%   - algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
%   - algoParams.lambda = 0.05; % Regularization parameter
%   - algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
%   - algoParams.TRY_PERIODIC_RESIDUAL = 0; % Take advantage of periodic residual if uniform TEs (will change range_fm)
%   - algoParams.residual: in case we pre-computed the fit residual (mostly for testing)
%
% Output: structure outParams
%   - outParams.species(ii).name: name of the species (taken from algoParams)
%   - outParams.species(ii).amps: estimated water/fat images, size [nx,ny,ncoils]
%   - outParams.r2starmap: R2* map (in s^{-1}, size [nx,ny])
%   - outParams.fieldmap: field map (in Hz, size [nx,ny])
%
%
% Author: Diego Hernando
% Date created: August 5, 2011
% Date last modified: November 10, 2011
% Date last modified: August 15, 2016 (BA Hoff)

function outParams = fw_i3cm1i_3pluspoint_hernando_graphcut(imData,algPar)

% Check validity of params, and set default algorithm parameters if not provided
[validParams,algPar] = checkParamsAndSetDefaults_graphcut(imData,algPar);
if validParams==0
    disp('Exiting -- data not processed');
end

% Get data dimensions
d = size(imData.images);
if (length(d)<5) || (d(5)<3)
    disp('Not enough TE images. Need at least 3.');
    outParams = [];
    return;
end
% [sx,sy,sz,C,N] = size(imData.images);

% If more than one channel, coil combine
if d(4) > 1
    disp('Multi-coil data: coil-combining');
    imData.images = coilCombine(imData.images);
end

% If precession is clockwise (positive fat frequency) simply conjugate data
if imData.PrecessionIsClockwise <= 0
    imData.images = conj(imData.images);
    imData.PrecessionIsClockwise = 1;
end

% Check spatial subsampling option (speedup ~ quadratic SUBSAMPLE parameter)
if algPar.SUBSAMPLE > 1
    timg = imData.images;
    START = round(algPar.SUBSAMPLE/2);
    [sx,sy] = size(timg(:,:,1,1,1));
    allX = 1:sx;
    allY = 1:sy;
    subX = START:algPar.SUBSAMPLE:sx;
    subY = START:algPar.SUBSAMPLE:sy;
    imData.images = timg(subX,subY,:,:,:);
end

% Let's get the residual. If it's not already in the params, compute it
if isfield(algPar,'residual')
    % Grab the residual from the params structure
    residual = algPar.residual;
else
    
    UNIFORM_TEs = 0;
    if (~isfield(algPar.TRY_PERIODIC_RESIDUAL) || algPar.TRY_PERIODIC_RESIDUAL)
        dTE = diff(imData.TE);
        UNIFORM_TEs = sum(abs(dTE - dTE(1))) < 1e-6;
    end
    
    % Compute the residual
    if UNIFORM_TEs == 1 % TEST DH* 090801
        % Find out the period in the residual (Assuming uniformly spaced samples)
        dt = imData.TE(2)-imData.TE(1);
        period = abs(1/dt);
        NUM_FMS_ORIG = algPar.NUM_FMS;
        range = diff(algPar.range_fm);
        algPar.NUM_FMS = ceil(algPar.NUM_FMS/range*period);
        algPar.range_fm = [0 period*(1-1/(algPar.NUM_FMS))];
        residual = computeResidual(imData,algPar);
        num_periods = ceil(range/period/2);
        algPar.NUM_FMS = 2*num_periods*algPar.NUM_FMS;
        residual = repmat(residual,[2*num_periods 1 1]);
        algPar.range_fm = [-num_periods*period (num_periods*period-period/NUM_FMS_ORIG)];
    else
        % If not uniformly spaced TEs, get the residual for the whole range
        residual = computeResidual(imData,algPar);
    end
    
end

% Setup the estimation, get the lambdamap,...
%   Spatially-varying regularization.  The LMAP_POWER applies to the
%       sqrt of the curvature of the residual, and LMAP_POWER=2 yields
%       approximately uniform resolution.
%   LMAP_EXTRA: Extra flexibility for including prior knowledge into
%       regularization. For instance, it can be used to add more smoothing
%       to noise regions (by adding, eg a constant LMAP_EXTRA), or even to
%       add spatially-varying smoothing as a function of distance to
%       isocenter...
fms = linspace(algPar.range_fm(1),algPar.range_fm(2),algPar.NUM_FMS);
dfm = fms(2)-fms(1);
lmap = getQuadraticApprox( residual, dfm );
lmap = (sqrt(lmap)).^algPar.LMAP_POWER;
lmap = lmap + mean(lmap(:))*algPar.LMAP_EXTRA;

% Initialize the field map indices
cur_ind = ceil(length(fms)/2)*ones(size(imData.images(:,:,1,1,1)));

% This is the core of the algorithm
fm = graphCutIterations(imData,algPar,residual,lmap,cur_ind);

% If we have subsampled (for speed), let's interpolate the field map
if algPar.SUBSAMPLE>1
    fmlowres = fm;
    [SUBX,SUBY] = meshgrid(subY(:),subX(:));
    [ALLX,ALLY] = meshgrid(allY(:),allX(:));
    fm = interp2(SUBX,SUBY,fmlowres,ALLX,ALLY,'*spline');
    lmap = interp2(SUBX,SUBY,lmap,ALLX,ALLY,'*spline');
    fm(isnan(fm)) = 0;
    imData.images = images0;
end
algPar.lmap = lmap;

% Now take the field map fm and get the rest of the estimates
r2starmap = zeros(d(1:3));
for ka=1:size(imData.images,6)
    
    curParams = imData;
    curParams.images = imData.images(:,:,:,:,:,ka);
    
    
    if algPar.range_r2star(2)>0
        % DH* 100422 use fine R2* discretization at this point
        algPar.NUM_R2STARS = round(algPar.range_r2star(2)/2)+1;
        r2starmap(:,:,ka) = estimateR2starGivenFieldmap( curParams, algPar, fm );
    else
        r2starmap(:,:,ka) = zeros(size(fm));
    end
    
    if algPar.DO_OT ~= 1
        % If no Optimization Transfer, just get the water/fat images
        amps = decomposeGivenFieldMapAndDampings( curParams,algPar, fm,r2starmap(:,:,ka),r2starmap(:,:,ka) );
        waterimage = squeeze(amps(:,:,1,:));
        fatimage = squeeze(amps(:,:,2,:));
        w(:,:,:,ka) = waterimage;
        f(:,:,:,ka) = fatimage;
    end
end

% If Optimization Transfer is requested, do it now
if algPar.DO_OT == 1
    imData.fmGC = fm;
    
    algPar.OT_ITERS = 10;
    algPar.fieldmap = fm;
    algPar.r2starmap = r2starmap;
    algPar.lambdamap = sqrt(algPar.lambda*lmap);
    
    outParams = fw_i3cm0i_3plusploint_hernando_optimtransfer( imData, algPar  );
    fm = outParams.fieldmap;
    
    
    % Now re-estimate the R2* map and water/fat images
    if algPar.range_r2star(2)>0
        algPar.NUM_R2STARS = round(algPar.range_r2star(2)/2)+1;
        r2starmap(:,:,ka) = estimateR2starGivenFieldmap( curParams,algPar, fm );
    else
        r2starmap(:,:,ka) = zeros(size(fm));
    end
    
    amps = decomposeGivenFieldMapAndDampings( curParams,algPar, fm,r2starmap(:,:,ka),r2starmap(:,:,ka) );
    waterimage = squeeze(amps(:,:,1,:));
    fatimage = squeeze(amps(:,:,2,:));
    w(:,:,:,ka) = waterimage;
    f(:,:,:,ka) = fatimage;
    
end

% Put results in outParams structure
try
    outParams.species(1).name = algPar.species(1).name;
    outParams.species(2).name = algPar.species(2).name;
catch
    outParams.species(1).name = 'water';
    outParams.species(2).name = 'fat';
end

outParams.species(1).amps = w;
outParams.species(2).amps = f;
outParams.r2starmap = r2starmap;
outParams.fieldmap = fm;




