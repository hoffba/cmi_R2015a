function [DT,Dtrace] = DT_All(dwi,diffdir,b,mask,fnout,fov)
% DT = DT_All(dwi,diffdir,b);
% Calculates diffusion tensor. Allows single or multiple b-values.
% Inputs:
%   dwi = [nx,ny,nz,ndwi] 4D image matrix
%   diffdir = [3 x ndwi] diffusion directions
%   b = b-values corresponding to dwi images (s/mm^2)
% (optional):
%   mask = binary mask applied to limit analysis to specified region
%   fnout = fulle file name for saving results
% Output:
%   DT = [nx,ny,nz,6] diffusion tensor images - Dxx,Dyy,Dzz,Dxy,Dxz,Dyz
%   Dtrace = 3D trace ADC
% This program calculates the Diffusion Tensor for Multiple B values and
% Multiple directions.

d = size(dwi);
if length(d)~=4
    error('Input image matrix must be 4D.');
elseif d(4)<7
    error('Not enough diffusion directions.');
end
if size(diffdir,1)~=3
    error('Invalid diffusion directions matrix');
elseif size(diffdir,2)~=d(4)
    error('Number of diffusion directions does not match DWI');
end
if length(b)~=d(4)
    error('Number of b-values does not match DWI.');
end
if (nargin>=5) && ~isempty(fnout) && ischar(fnout)
    fdir = fileparts(fnout);
    if ~exist(fdir,'dir')
        fnout = '';
    end
end
if (nargin<6) || ~isnumeric(fov) || (length(fov)~=3)
    fov = d(1:3);
end

% Find b0 image:
ind = b==0;
b0 = mean(dwi(:,:,:,ind),4);
dwi(:,:,:,ind) = [];
b(ind) = [];
diffdir(:,ind) = [];

% Initialize the b-matrix:
n = length(b);
B = zeros(n,6);
for i = 1:n
    tB = b(i)*diffdir(:,i)*diffdir(:,i)';
    B(i,:) = [diag(tB)',2*tB([2,3,6])];
end
Bi = pinv(B);

% Prep the images:
for i = 1:n
    dwi(:,:,:,i) = -log(dwi(:,:,:,i)./b0);
end

% Solve diffusion equation for each voxel:
DT = zeros([d(1:3),6]);
n = prod(d(1:3));
mask = find(mask);
[i,j,k] = ind2sub(d(1:3),mask);
hw = waitbar(0);
for ii = 1:length(mask)
    waitbar(mask(ii)/n,hw,sprintf('Calculating Diffusion Tensor: %u / %u',mask(ii),n));
    DT(i(ii),j(ii),k(ii),:) = Bi*squeeze(dwi(i(ii),j(ii),k(ii),:));
end
delete(hw);

% Remove erroneous voxels:
ind = repmat(any(isinf(DT) | isnan(DT),4),1,1,1,6);
% Diagonal D values cannot be negative:
ind = ind | cat(4,repmat(any(DT(:,:,:,1:3)<0,4),1,1,1,3),zeros([d(1:3),3]));
DT(ind) = 0;

if nargout>1
    Dtrace = sum(DT(:,:,:,1:3),4)/3;
end

if ~isempty(fnout)
    saveMHD(fnout,cat(4,DT,Dtrace),{'Dxx','Dyy','Dzz','Dxy','Dxz','Dyz','Dtrace'},fov)
end