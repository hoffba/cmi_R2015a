% Function: computeResidual
%
% Description: compute fit residual for water/fat imaging
% 
% Parameters:
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,nz,nTE]
%%   - imDataParams.TEs: echo times (in seconds)
%%   - imDataParams.fieldStrength: (in Tesla)
%%
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
%
% Returns: 
%  - residual: the residual, of size NUM_FMS X sx X sy
%
% Author: Diego Hernando
% Date created: August 13, 2011
% Date last modified: December 8, 2011

function residual = computeResidual3d( imDataParams, algoParams )

images = imDataParams.images;

try
  precessionIsClockwise = imDataParams.PrecessionIsClockwise;
catch
  precessionIsClockwise = 1;
end

  
% If precession is clockwise (positive fat frequency) simply conjugate data
if precessionIsClockwise <= 0 
  imDataParams.images = conj(imDataParams.images);
  imDataParams.PrecessionIsClockwise = 1;
end

gyro = 42.58;
deltaF = [0 ; gyro*(algoParams.species(2).frequency(:) - algoParams.species(1).frequency(1))*(imDataParams.FieldStrength)];
relAmps = algoParams.species(2).relAmps;
range_fm = algoParams.range_fm;
t = imDataParams.TE;
NUM_FMS = algoParams.NUM_FMS; 
range_r2star = algoParams.range_r2star;
NUM_R2STARS = algoParams.NUM_R2STARS;
r2s =linspace(range_r2star(1),range_r2star(2),NUM_R2STARS);
psis = linspace( range_fm(1),range_fm(2),NUM_FMS );

[sx,sy,sz,N] = size(images);
np = sx*sy*sz;
images = reshape(images,np,N).';

% Get masked values:
if isfield(imDataParams,'mask')
    mask = imDataParams.mask;
else
    mask = true(sx,sy,sz);
end

% Get VARPRO-formulation matrices for given echo times and chemical shifts 
Phi = getPhiMatrixMultipeak(deltaF,relAmps,t);
residual = zeros(NUM_FMS,np);
for k = 1:NUM_FMS
    fprintf('FMS: %u/%u\n',k,NUM_FMS);
    for kr = 1:NUM_R2STARS
        % Compute projector matrix:
        Psi = diag(exp((1i*2*pi*psis(k) - r2s(kr))*t));
        P = (eye(N)-Psi*Phi*pinv(Psi*Phi));
        % Sum over echoes and coils:
        tR = sum(abs(P*images(:,mask)).^2,1);
        % Find min over R2star
        if kr==1
            R = tR;
        else
            R = min(R,tR);
        end
    end
    residual(k,mask) = R;
end
residual = reshape(residual,[NUM_FMS,sx,sy,sz]);



