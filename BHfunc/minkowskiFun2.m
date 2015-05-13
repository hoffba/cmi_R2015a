%
function [MF,labels] = minkowskiFun2(BW,n,mask)
% Calculates regional minkowski functionals using a moving window method
% Inputs:
%   BW      = 3D binary matrix
%   n       = moving window radius
%   mask    = binary image mask for VOI
% Output: fimg = 4D minkowski functional maps 
%                   OR [#thresh x nMF] matrix of results for global values
%           * 2D : Area, Perimeter, Euler-Poincare
%           * 3D : Volume, Surface Area, Mean Breadth, Euler-Poincare
%         labels = cell array of MF names

d = size(BW);
if (nargin<3)
    mask = [];
elseif ~((ndims(BW)==ndims(mask)) && all(d==size(mask)))
    warning('Mask dimensions do not match, proceeding without a mask ...');
    mask = [];
end
if nargin<2
    R = [4,4,0];
else
    R = zeros(1,3);
    R(1:length(n)) = n;
end
nD = nnz(R);
switch nD
    case 2
        func = @calcMF2D;
        labels = {'Area','Perimeter','Euler'};
        nmf = 3;
    case 3
        func = @calcMF3D;
        labels = {'Volume','SurfaceArea','Mean Breadth','Euler'};
        nmf = 4;
    otherwise
        error(['Invalid window size: [',num2str(n),']']);
end

if any(isinf(n))
    % Perform analysis on entire image
    if ~isempty(mask)
        BW = BW & mask;
    end
    MF = feval(func,BW);
else
    % Perform analysis on moving window
    if isempty(mask)
        ind = 1:prod(d(1:3));
    else
        ind = find(mask);
    end
    ntot = length(ind);
    disp(['Number of iterations: ',num2str(ntot)])

    MF = zeros([d,nmf]);
    hw = waitbar(0,'Calculating local Minkowski Functionals ...');
    t = tic;
    for i = 1:ntot
        if mod(i,10000)==0
            disp(['  ',num2str(ntot-i),' left (',...
                  datestr(toc(t)*(ntot-i)/(i*60*60*24),'DD:HH:MM:SS'),')'])
        end
        [ii,jj,kk] = ind2sub(d(1:3),ind(i));
        iBW = BW(max(1,ii-n(1)):min(d(1),ii+n(1)),...
                 max(1,jj-n(2)):min(d(2),jj+n(2)),...
                 max(1,kk-n(3)):min(d(3),kk+n(3)));
        MF(ii,jj,kk,:) = feval(func,iBW);
        waitbar(i/ntot,hw,['Complete: ',num2str(i),'/',num2str(ntot)]);
    end
    delete(hw);
    toc(t);
end

% Raw 2D MF calculation
function p = calcMF2D(BW)
p = [ imAreaDensity(BW) , ...
      imPerimeterDensity(BW) ,...
      imEuler2dDensity(BW) ];
  
% Raw 3D MF calculation
function p = calcMF3D(BW)
p = [ imVolumeDensity(BW) ,...
      imSurfaceDensity(BW) ,...
      imMeanBreadth(BW) ,...
      imEuler3dDensity(BW)];



