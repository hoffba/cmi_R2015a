% Function: computeResidual
%
% Description: compute fit residual for water/fat imaging
% 
% Parameters:
% Input: structures imData and algPar
%   - imData.images: acquired images, array of size[nx,ny,nz,ncoils,nTE]
%   - imData.TEs: echo times (in seconds)
%   - imData.fieldStrength: (in Tesla)
%
%   - algPar.species(ii).frequency = frequency shift in ppm of each peak within species ii
%   - algPar.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%   Example
%      - algPar.species(1).name = 'water' % Water
%      - algPar.species(1).frequency = [0] 
%      - algPar.species(1).relAmps = [1]   
%      - algPar.species(2).name = 'fat' % Fat
%      - algPar.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%      - algPar.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
% 
%   - algPar.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
%   - algPar.range_r2star = [0 0]; % Range of R2* values
%   - algPar.NUM_R2STARS = 1; % Numbre of R2* values for quantization
%   - algPar.range_fm = [-400 400]; % Range of field map values
%   - algPar.NUM_FMS = 301; % Number of field map values to discretize
%   - algPar.NUM_ITERS = 40; % Number of graph cut iterations
%   - algPar.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
%   - algPar.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
%   - algPar.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
%   - algPar.lambda = 0.05; % Regularization parameter
%   - algPar.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
%   - algPar.TRY_PERIODIC_RESIDUAL = 0; % Take advantage of periodic residual if uniform TEs (will change range_fm)  
%   - algPar.residual: in case we pre-computed the fit residual (mostly for testing) 
%
% Returns: 
%  - residual: the residual, of size NUM_FMS X sx X sy
%
% Author: Diego Hernando
% Date created: August 13, 2011
% Date last modified: December 8, 2011

function residual = computeResidual(imData,algPar)

% Size should be [ sx,sy,sz, C coils , N echoes , a acquisitions ]
[d(1),d(2),d(3),d(4),d(5),d(6)] = size(imData.images);

% If precession is clockwise (positive fat frequency) simply conjugate data
if isfield(imData,'PrecessionIsClockwise') && (imData.PrecessionIsClockwise <= 0)
    imData.images = conj(imData.images);
    imData.PrecessionIsClockwise = 1;
end

gyro = 42.58;
deltaF = [0 ; gyro*(algPar.species(2).frequency(:) - algPar.species(1).frequency(1))*(imData.FieldStrength)];
relAmps = algPar.species(2).relAmps;
range_fm = algPar.range_fm;

% Get VARPRO-formulation matrices for given echo times and chemical shifts
Phi = getPhiMatrixMultipeak(deltaF,relAmps,imData.TE);
%Phi(:,1) = 0;
% iPhi = pinv(Phi'*Phi);
% A = Phi*iPhi*Phi';

psis = linspace( range_fm(1),range_fm(2),algPar.NUM_FMS );

% Compute residual
r2s =linspace(algPar.range_r2star(1),algPar.range_r2star(2),algPar.NUM_R2STARS);

% Precompute all projector matrices (one per field value) for VARPRO
P = [];
for kr = 1:algPar.NUM_R2STARS
    P1 = [];
    for k = 1:algPar.NUM_FMS
        Psi = diag(exp(1j*2*pi*psis(k)*imData.TE - abs(imData.TE)*r2s(kr)));
        P1 = [P1;(eye(d(2))-Psi*Phi*pinv(Psi*Phi))];
%         P1 = [P1;(eye(d(5))-Psi*Phi*pinv(Psi*Phi))];
    end
    P(:,:,kr) = P1;
end

% Compute residual for all voxels and all field values
% Note: the residual is computed in a vectorized way, for increased speed
residual = zeros(algPar.NUM_FMS,d(1),d(2));
% $$$   r2array = zeros(NUM_FMS,d(1),d(2),num_acq);

% Go line-by-line in the image to avoid using too much memory, while
% still reducing the loops significantly
for ka = 1:d(6) % # of acquisitions
    for ky=1:d(2)
        temp = reshape(squeeze(permute(imData.images(:,ky,:,:,:,ka),[1 2 3 5 4])),[d(1) d(2)*d(2)]).';
        temp = reshape(temp,[d(2) d(1)*d(2)]);
        temp3 = zeros(d(2),algPar.NUM_R2STARS);
        for kr = 1:algPar.NUM_R2STARS
            temp2 = reshape(sum(abs(reshape(P(:,:,kr)*temp,[d(2) d(2)*algPar.NUM_FMS*d(1)])).^2,1),[algPar.NUM_FMS d(2)*d(1)]).';
            temp3(:,kr) = sum(reshape(temp2(:,:,kr),[d(2) algPar.NUM_FMS*d(1)]),1);
        end
        [mint3,~] = min(temp3,[],2);
        
        residual(:,:,ky) = squeeze(squeeze(residual(:,:,ky)).' + reshape(mint3,[d(1) algPar.NUM_FMS])).';
        % $$$       r2array(:,:,ky,ka) = (reshape(r2s(imint3),[sx NUM_FMS])).';
    end
end

