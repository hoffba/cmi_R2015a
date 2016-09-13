function [DT1,Dtrc] = DT_All_mk1(dwi,diffdir,b,mask)
% DT1 = DT_All_mk1(dwi,diffdir,b);
% Calculates diffusion tensor. Allows single or multiple b-values.
% Inputs:
%   dwi = [nx,ny,nz,ndwi] 4D image matrix
%   diffdir = [3 x ndwi] diffusion directions
%   (NOTE: assumes that "dwi" does NOT contain "trace" images, but includes b=0)
%   b = b-values corresponding to dwi images (s/mm^2)
% (optional):
%   mask = binary mask applied to limit analysis to specified region
%   fnout = fulle file name for saving results
% Output:
%   DT1 = [nx,ny,nz,6] diffusion tensor images - Dxx,Dyy,Dzz,Dxy,Dxz,Dyz
%   Dtrace = 3D trace ADC
% This program calculates the Diffusion Tensor for Multiple B values and
% Multiple directions.
%
% original code source: "DT_ALL" script by Ben Hoff (July 2016)
% MD 20160728: vectorized DT1 calculation
% MD 20160729: attempted to vectorize "filtering" of unphysical voxels, 
% but got different results from Ben's 3D filter for "diagonal"<0 part??
%

szd = size(dwi);
if length(szd)~=4
    error('Input image matrix must be 4D.');
elseif szd(4)<7
    error('Not enough diffusion directions.');
end
if size(diffdir,1)~=3
    error('Invalid diffusion directions matrix');
elseif size(diffdir,2)~=szd(4)
    error('Number of diffusion directions does not match DWI');
end
if length(b)~=szd(4)
    error('Number of b-values does not match DWI.');
end

% Find b0 image:
ind = b==0;
dwib0 = mean(dwi(:,:,:,ind),4);
dwi(:,:,:,ind) = [];
b(ind) = [];
diffdir(:,ind) = [];

% Initialize the b-matrix:
nd = length(b); %debug % Ndif-dirs x Nb >0
B = zeros(nd,6);

for i = 1:nd % for each diffusion direction and b>0
    tB = b(i)*diffdir(:,i)*diffdir(:,i)'; %assuming "normalized" diff-n vectors
    B(i,:) = [diag(tB)',2*tB([2,3,6])];
end
Bi = pinv(B); % size of B'(transpose)
%size(Bi)  %debug

% Solve diffusion equation for each voxel:
DT1 = zeros([szd(1:3),6]);
n = prod(szd(1:3));
% size(DT) %debug

% MD 20160728: vectorized DT1 calculation
DTv = reshape(DT1, [n 6])';
dwiv = reshape(dwi, [n nd])';

imsk1d = find(mask(:)>0); %1D mask index
% size(repmat(dwib0(imsk1d)', nd,1)) %debug
% size(dwiv(:,imsk1d))% debug

dwiv(:,imsk1d) = -log(dwiv(:,imsk1d)./(repmat(dwib0(imsk1d)', nd,1)+eps));
% dwiv = -log(dwiv./(repmat(dwib0(:)', nd,1)+eps));
dwiv(:,mask(:)<1)=0; % set hard zeros for the mask
DTv(:,imsk1d) = Bi*dwiv(:,imsk1d); % MD 20160728: vectorized DT1 calculation

% size(DTv) % debug
% size(dwiv) % debug
% size(imsk1d) % debug

% % Filter unphysical voxels in 2D -- not same as in 3D below??
%  DTvd = DTv(1:3,:);
%  % size(DTvd) %debug
%  DTvd(DTvd(:)<0)=0; DTv(1:3,:)=DTvd; % DTii<0 not allowed
%  clear DTvd;
%  DTv(isnan(DTv(:))>0 | isinf(DTv(:))>0)=0;

DT1 = reshape(DTv',[szd(1:3) 6]); % make a tensor
clear DTv;
% size(DT1) % debug
% sum(DT1(:)) %debug

% % Remove erroneous voxels in 3D:
ind1 = repmat(any((isinf(DT1) | isnan(DT1)),4),[1 1 1 6]);
% Diagonal D values cannot be negative:
ind1 = ind1 | cat(4,repmat(any(DT1(:,:,:,1:3)<0,4),[1 1 1 3]),zeros([szd(1:3),3]));
DT1(ind1) = 0;

%sum(DT1(:)) %debug

if nargout>1
    Dtrc = sum(DT1(:,:,:,1:3),4)/3;
end

return
