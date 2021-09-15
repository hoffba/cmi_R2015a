function [Vmax,Nmax] = MTHT3D(im,s,no,beta,alpha,c,alfa)
%%
%	INPUT:
%       im - 3D input image (image normalized to between 0,1)
%       o - Orientation angles
%		s - Scale
%		no - Number of orientations
%		c - Parameter for vesselness
% 		beta - Parameter for vesselness
% 		alfa - parameter for the neuriteness
%
%   OUTPUT:
%       Vmax - Vesselness Enhanced Image
%		Nmax - Neuriteness Enhanced Image
%
%   USAGE:
%       [Vmax,Nmax] = MTHT2D(im,o,s,no,beta,c,alfa)
%
%   AUTHOR:
%       Sundaresh Ram
% 		Galban Lab
% 		Dept. Of Radiology & Dept. of BME
% 		University of Michigan
% 		e-mail : sundarer@umich.edu
%
%   VERSION:
%       0.1 - 02/21/2020 First implementation
%
%
Nmax = 0;
TH = zeros(size(im,1),size(im,2),size(im,3),no);
V = zeros(size(im,1),size(im,2),size(im,3),length(s));
%N = zeros(size(im,1),size(im,2),size(im,3),length(s));
%% Point distribution on the sphere of unit radius
[p,ps,~] = SurfacesSpiralPoints3D(no);
orients = ps;
for i = 1:length(s)
    fprintf(['Scale: ' num2str(s(i)),'\n']);
    for j = 1:length(p)
        %fprintf(['Scale: ' num2str(s(i)), '  orientation: ' num2str(ps(j)) '\n']);
        se = Line3D(s(i),p(j,:),1,1,1); %3D line structuring element 
        TH(:,:,:,j) = im - imopen(im,se); 
    end
    %figure, imshow3D(TH);
    %% Tensor
    T = TensorForOrientedQuadratureFilters3D(TH,orients); 
    %% Eigen Matrix - values and vectors
    [L1,L2,L3] = EigenMatrix3x3M(T);
    [L1,L2,L3] = EigenSort3x3M(L1,L2,L3); % sort
    %% Vesselness
    imf1 = Vesselness3D(im,beta,alpha,c,L1,L2,L3);
    %% Neuritenees
    %imf2 = Neuritenees3D(alfa,L1,L2,L3);
    %% Table of all enhanced image in each scale
    V(:,:,:,i) = imf1;
    %N(:,:,:,i) = imf2; 
end
%% Calculate Maximum Image Over the Scales
if length(s)>1
   Vmax = squeeze(max(V,[],4));
  % Nmax = squeeze(max(N,[],4));
else
   Vmax = squeeze(V);
  % Nmax = squeeze(N);
end
end %End of function