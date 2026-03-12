function Vmax = MTHT3D_BH(im,s,no,beta,alpha,c)
%%
%	INPUT:
%       im - 3D input image (image normalized to between 0,1)
%       o - Orientation angles
%		s - Scale
%		no - Number of orientations
%		c - Parameter for vesselness
% 		beta - Parameter for vesselness
%
%   OUTPUT:
%       Vmax - Vesselness Enhanced Image
%
%   USAGE:
%       [Vmax,Nmax] = MTHT2D(im,o,s,no,beta,c)
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

d = size(im,1:3);
TH = single(zeros([d,no]));
Vmax = zeros(d);
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
    [L1,L2,L3] = EigenMatrix3x3M(T); clear T
    [L1,L2,L3] = EigenSort3x3M(L1,L2,L3); % sort
    %% Vesselness
    Vmax = max(Vesselness3D(im,beta,alpha,c,L1,L2,L3),Vmax);
end