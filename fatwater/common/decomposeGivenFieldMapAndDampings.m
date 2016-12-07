% Function: decomposeGivenFieldMapAndDampings
%
% Description: estimate water/fat images given the nonlinear parameters
% 
% Parameters:
%% Input: structures imDataParams and algoParams
%%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
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
%
%  - fieldmap: the estimated B0 field map
%  - r2starWater: the estimated water R2* map
%  - r2starFat: the estimated fat R2* map
%
% Returns: 
%  - amps: the amplitudes for all chemical species and coils
%  - remerror: fit error norm
%
% Author: Diego Hernando
% Date created: 
% Date last modified: August 18, 2011

function [amps,remerror] = decomposeGivenFieldMapAndDampings( imDataParams,algoParams,fieldmap,r2starWater,r2starFat )

gyro = 42.58;

try
    ampW = algoParams.species(1).relAmps;
catch
    fprintf('Setting ampW = 1 (default)');
    ampW = 1.0;
end

% If precession is clockwise (positive fat frequency) simply conjugate data
if ~isfield(imDataParams,'precessionIsClockwise') || (imDataParams.PrecessionIsClockwise <= 0)
    imDataParams.images = conj(imDataParams.images);
    imDataParams.PrecessionIsClockwise = 1;
end

deltaF = [0 ; gyro*(algoParams.species(2).frequency(:) - algoParams.species(1).frequency(1))*(imDataParams.FieldStrength)];
relAmps = algoParams.species(2).relAmps;
images = imDataParams.images;
t = imDataParams.TE;

[sx,sy,sz,C,N,~] = size(images);
images = reshape(imDataParams.images,[],N).';

% Apply only to masked indices:
if isfield(imDataParams,'mask')
    mchk = true;
    fieldmap = fieldmap(imDataParams.mask);
    images = images(:,imDataParams.mask);
else
    mchk = false;
    fieldmap = fieldmap(:);
end
np = length(fieldmap); % Number of spatial points being calculated

relAmps = reshape(relAmps,1,[]);

B1 = zeros(N,2);
B = zeros(N,2);
for n=1:N
    B1(n,:) = [ampW*exp(1i*2*pi*deltaF(1)*t(n)),sum(relAmps(:).*exp(1i*2*pi*deltaF(2:end)*t(n)))];
end

A = B\images;

if mchk
    amps = zeros();
    amps() = amps
else
    amps = reshape(A,sx,sy,sz,C,2);
end



% !! OLD !!

remerror = zeros(sx,sy,sz);
amps = nan(sx,sy,sz,2);
for kx =1:sx
    for ky=1:sy
        for kz=1:sz
            s = reshape( squeeze(images(kx,ky,kz,:,:)), [C N]).';

            B(:,1) = B1(:,1).*exp(1i*2*pi*fieldmap(kx,ky,kz)*t(:) - r2starWater(kx,ky,kz)*t(:));
            B(:,2) = B1(:,2).*exp(1i*2*pi*fieldmap(kx,ky,kz)*t(:) - r2starFat(kx,ky,kz)*t(:));

            amps(kx,ky,kz,:) = B\s;

            if nargout == 2
                remerror(kx,ky) = norm(s - B*squeeze(amps(kx,ky,:,:)),'fro');
            end
        end
    end
end


