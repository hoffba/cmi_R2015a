% Function: estimateR2starGivenFieldmap
%
% Description: estimate R2* map, given the fieldmap
%
% Parameters:
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%%   - imDataParams.TEs: echo times (in seconds)
%%   - imDataParams.fieldStrength: (in Tesla)
%%
%%   - algoParams.species(ii).name = string containing the name of the species
%%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%%   Example
%%      - algoParams.species(1).name = 'water' % Water
%%      - algoParams.species(1).frequency = [0]
%%      - algoParams.species(1).relAmps = [1]
%%      - algoParams.species(2).name = 'fat' % Fat
%%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
%%   - algoParams.range_r2star = [0 0]; % Range of R2* values
%%   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
%%
%  - fm: the estimated B0 field map
%
% Returns:
%  - r2starmap: the estimated R2* map
%  - residual: fit error residual
%
% Author: Diego Hernando
% Date created: 2009
% Date last modified: August 18, 2011


function [r2starmap,residual] = estimateR2starGivenFieldmap3d ( imDataParams, algoParams, fm )

if isfield(imDataParams,'PrecessionIsClockwise')
    precessionIsClockwise = imDataParams.PrecessionIsClockwise;
else
    precessionIsClockwise = 1;
end
% If precession is clockwise (positive fat frequency) simply conjugate data
if precessionIsClockwise <= 0
    imDataParams.images = conj(imDataParams.images);
    imDataParams.PrecessionIsClockwise = 1;
end

[sx,sy,sz,C,N,num_acqs] = size(imDataParams.images);
images = permute(reshape(imDataParams.images,[sx*sy*sz,C,N,num_acqs]),[3,1,2,4]);

if isfield(imDataParams,'mask')
    mchk = true;
    fm = fm(imDataParams.mask);
    images = images(:,imDataParams.mask,:,:);
else
    mchk = false;
    fm = fm(:);
end
np = length(fm); % Number of spatial points being calculated

range_r2star = algoParams.range_r2star;
NUM_R2STARS = algoParams.NUM_R2STARS;
gyro = 42.58;
deltaF = [0 ; gyro*(algoParams.species(2).frequency(:) - algoParams.species(1).frequency(1))*(imDataParams.FieldStrength)];
relAmps = algoParams.species(2).relAmps;
t = imDataParams.TE(:);

% Undo effect of field map
fmcorr = nan(N,np);
for kt = 1:N
    fmcorr(kt,:) = exp(-1i*2*pi*fm*t(kt));
end

r2s = linspace(range_r2star(1),range_r2star(2),NUM_R2STARS);
Phi = getPhiMatrixMultipeak(deltaF,relAmps,t);
R = inf(np,1); iminres = zeros(np,1);
for kr = 1:NUM_R2STARS
    % Compute projector matrix:
    Psi = diag(exp(-abs(t)*r2s(kr)));
    P = (eye(N)-Psi*Phi*pinv(Psi*Phi));
    tR = zeros(np,1);
    for ka = 1:num_acqs
        for kC = 1:C
            S = reshape(images(:,:,kC,ka),N,[]) .* fmcorr;
            % Sum over echoes, coils, and acquisitions:
            tR = tR + sum(abs(P*S).^2,1)';
        end
    end
    % Find min over R2star
    iR = tR<R;
    R(iR) = tR(iR);
    iminres(iR) = kr;
end
if mchk
    residual = zeros(sx,sy,sz);
    residual(imDataParams.mask) = R;
    r2starmap = zeros(sx,sy,sz);
    r2starmap(imDataParams.mask) = r2s(iminres);
else
    residual = reshape(R,[sx,sy,sz]);
    r2starmap = r2s(iminres);
end

