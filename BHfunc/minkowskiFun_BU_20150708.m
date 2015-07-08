%
function [MF,labels] = minkowskiFun(BW,n,ind,voxsz)
% Calculates regional minkowski functionals using a moving window method
% Inputs:
%   BW      = 2D/3D binary matrix
%   n       = moving window radius
%   ind     =  (3D) indices of matrix locations to analyze
% Output:
%   MF     = [#pts x nMF] Minkowski Functionals 
%               * 2D : Area, Perimeter, Euler-Poincare
%               * 3D : Volume, Surface Area, Mean Breadth, Euler-Poincare
%   labels = cell array of MF names

BW = logical(BW);
d = size(BW);
if (nargin<3)
    ind = [];
elseif any((ind<1) | (ind>prod(d)))
    warning('Mask dimensions do not match, proceeding without a mask ...');
    ind = [];
end
if nargin<2
    n = inf([1,length(d)]);
end
gchk = any(isinf(n));
if gchk
    nD = length(d);
else
    nD = nnz(n);
end
if (nargin<4)
    voxsz = ones(nD);
else
    voxsz = voxsz(1:nD);
end

switch nD
    case 2
        if gchk
            func = @calcMF2D;
        else
            func = @calcMF2Drho;
        end
        labels = {'Area','Perimeter','Euler'};
        nmf = 3;
    case 3
        if gchk
            func = @calcMF3D;
        else
            func = @calcMF3Drho;
        end
        labels = {'Volume','SurfaceArea','Mean Breadth','Euler'};
        nmf = 4;
    otherwise
        error(['Invalid window radius: [',num2str(n),']']);
end

if gchk
    % Perform global analysis
    if nD == length(d)
        MF = feval(func,BW,voxsz);
    else
        warning('Dimensionality does not match.');
        MF = [];
    end
else
    % Perform analysis on moving window
    if isempty(ind)
        ind = 1:prod(d);
    end
    ntot = length(ind);
    disp(['Number of iterations: ',num2str(ntot)])

    MF = zeros([ntot,nmf]);
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
        MF(i,:) = feval(func,iBW,voxsz);
        waitbar(i/ntot,hw,['Complete: ',num2str(i),'/',num2str(ntot)]);
    end
    delete(hw);
    toc(t);
end

% Raw 2D MF calculation
function p = calcMF2D(BW,voxsz)
p = [ imAreaEstimate(BW,voxsz) , ...
      imPerimeterEstimate(BW,voxsz) ,...
      imEuler2dEstimate(BW) ];
  
% Raw 3D MF calculation
function p = calcMF3D(BW,voxsz)
p = [ imVolumeEstimate(BW,voxsz) ,...
      imSurfaceEstimate(BW,voxsz) ,...
      imMeanBreadth(BW,voxsz) ,...
      imEuler3dEstimate(BW)];
  
% Raw 2D MF calculation
function p = calcMF2Drho(BW,voxsz)
p = [ imAreaDensity(BW,voxsz) , ...
      imPerimeterDensity(BW,voxsz) ,...
      imEuler2dDensity(BW) ];
  
% Raw 3D MF calculation
function p = calcMF3Drho(BW,voxsz)
p = [ imVolumeDensity(BW,voxsz) ,...
      imSurfaceDensity(BW,voxsz) ,...
      imMeanBreadth(BW,voxsz) ,...
      imEuler3dDensity(BW)];



