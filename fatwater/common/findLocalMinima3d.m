% Function name: findLocalMinima
%
% Description: Find the local minima of the VARPRO residual at each
% voxel, and return them in a 3D array.
%
% Input arguments:
%   - residual: the 3D residual array (of size NUM_FMS x SX x SY)
%   - threshold: a threshold for the signal at each voxel to be even considered 
%     (voxels with signal lower than the threshold will not be assigned any local 
%     minima, as these will be meaningless anyway)
%   - masksignal: alternatively, one can provide a binary map specifying which 
%     voxels have non-negligible signal
%
% Output arguments:
%   - masksignal: the masksignal (same as input, or the one computed inside 
%     function if not provided)
%   - resLocalMinima: the main output of this function: a 3D array giving, at 
%     each voxel, the indices where the local minima of the residual are located.
%   - numMinimaPerVoxel: SXxSY map specifying how many local minima were found 
%     in the residual at each voxel
%
%Author: Diego Hernando
%Date created: March 18, 2008
%Date last update: March 18, 2008
%--------------------------------------------------------------------------

function [masksignal,resLocalMinima,numMinimaPerVoxel] = findLocalMinima3d( residual, threshold, masksignal )


% Get the dimensions
[~,sx,sy,sz] = size(residual);
% L = size(residual,1);
% sx = size(residual,2);
% sy = size(residual,3);

% Take finite differences in residual
dres = diff(residual,1,1);

maxres = max(residual,[],1);
minres = min(residual,[],1);

% Find voxels with signal
if nargin < 3
  sumres = sqrt(squeeze(sum(residual,1)));
  sumres = sumres/max(sumres(:));
  masksignal = sumres>threshold;
end

resLocalMinima = zeros(1,sx,sy,sz);
numMinimaPerVoxel = zeros(sx,sy,sz);
for i = find(masksignal)'
    [xx,yy,zz] = ind2sub([sx,sy,sz],i);
    tmx = maxres(1,xx,yy,zz);
    tmn = minres(1,xx,yy,zz);
    tmp = [0;squeeze(dres(:,xx,yy,zz))];
    tmp = find(tmp<0 & circshift(tmp,-1)>0 & residual(:,xx,yy,zz)<tmn+0.3*(tmx-tmn));
    numMinimaPerVoxel(xx,yy,zz) = numel(tmp);
    resLocalMinima(1:numMinimaPerVoxel(xx,yy,zz),xx,yy,zz) = tmp;
end

% % Loop through all voxels (this is MATLAB-inefficient and should be done in vectorized way)
% resLocalMinima = zeros(1,sx,sy,sz);
% numMinimaPerVoxel = zeros(sx,sy,sz);
% for kx=1:sx
%   for ky=1:sy
%     if masksignal(kx,ky) > 0
%       minres = min(residual(:,kx,ky));
%       maxres = max(residual(:,kx,ky));
% 
%       temp = [0;squeeze(dres(:,kx,ky))];
%       temp = temp<0 & circshift(temp,-1)>0 & residual(:,kx,ky)<minres+0.3*(maxres-minres);
%       
%       resLocalMinima(1:sum(temp),kx,ky) = find(temp);
%       numMinimaPerVoxel(kx,ky) = sum(temp);
%     end
%   end
% end

% $$$ x = 1:sx;
% $$$ y = 1:sy;
% $$$ [Y,X] = meshgrid(y,x);
% $$$ 
% $$$ for kx=1:sx
% $$$   for ky=1:sy
% $$$     if masksignal(kx,ky) == 0
% $$$ 
% $$$       noise_vals = 5:10:size(residual,1);
% $$$ 
% $$$       resLocalMinima(1:length(noise_vals),kx,ky) = noise_vals';
% $$$       numMinimaPerVoxel(kx,ky) = length(noise_vals);
% $$$     end
% $$$   end
% $$$ end
% $$$ 
