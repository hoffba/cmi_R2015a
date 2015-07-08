%
function [MF,labels] = minkowskiFun2(img,thresh,tmode,varargin)
% Calculates regional minkowski functionals using a moving window method
% Inputs:
%   img     = 2D/3D binary matrix
%   thresh  = image thresholds for binarization
%   tmode   = (p indicates percentile) <,>,==,etc.
%   Optional Name/Value pairs:
%       n       = moving window radius
%       ind     = (3D) indices of matrix locations to analyze
%       voxsz   = voxel size in 3D
%       mask    = binary image mask
%       defVal  = (logical) value assigned to voxels outside the mask
%       prog    = logical check to display waitbar/progress
% Output:
%   MF     = [#pts x nMF] Minkowski Functionals 
%               * 2D : Area, Perimeter, Euler-Poincare
%               * 3D : Volume, Surface Area, Mean Breadth, Euler-Poincare
%   labels = cell array of MF names

p = inputParser;
addRequired(p,'img',@isnumeric);
addRequired(p,'thresh',@isvector);
addRequired(p,'tmode',@ischar);
addParameter(p,'n',inf(1,3),@isvector);
addParameter(p,'ind',[],@isvector);
addParameter(p,'voxsz',ones(1,3),@isvector);
addParameter(p,'mask',[],@islogical);
addParameter(p,'defVal',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'prog',false,@(x)islogical(x)&&isscalar(x));
parse(p,img,thresh,tmode,varargin{:});
pp = p.Results;
d = size(pp.img);
if any((pp.ind<1) | (pp.ind>prod(d)))
    warning('Invalid matrix index ... proceeding with full image.');
    pp.ind = [];
end
gchk = any(isinf(pp.n));
if gchk
    nD = length(d);
else
    nD = nnz(pp.n);
end
if strncmp(pp.tmode,'p',1)
    pchk = true;
    pp.tmode(1) = [];
else
    pchk = false;
end
if ~ismember(pp.tmode,{'<','<=','>','>=','==','~='})
    error(['Threshold function (',pp.tmode,') is not valid.']);
end
if ~isempty(pp.mask)
    if ~all(size(pp.mask)==d)
        error('Mask dimensions do not fit image');
    end
    mchk = true;
else
    pp.mask = true(d);
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
        error(['Invalid window radius: [',num2str(pp.n),']']);
end
if pchk
    labels = [labels,{'Threshold'}];
end
nth = length(pp.thresh);

if gchk % Perform global analysis
    
    MF = nan(nth,nmf);
    
    % Convert percentiles to threshold values:
    if pchk
        pp.thresh = prctile(pp.img(pp.mask),pp.thresh);
    end
    if pp.prog, hw = waitbar(0,'Processing ...'); end
    for ith = 1:nth
        BW = eval(['pp.img ',pp.tmode, 'pp.thresh(ith);']);
        BW(~pp.mask) = pp.defVal;
        MF(ith,:) = feval(func,BW,pp.voxsz);
        if pp.prog, waitbar(ith/nth,hw); end
    end
    if pp.prog, delete(hw); end
    if pchk
        MF = [MF,pp.thresh(:)];
    end
    
else
    % Perform analysis on moving window
    if isempty(pp.ind)
        pp.ind = 1:prod(d);
    end
    ntot = length(pp.ind);
    disp(['Number of iterations: ',num2str(ntot)])

    MF = zeros(nth,nmf+pchk,ntot);
    hw = waitbar(0,'Calculating local Minkowski Functionals ...');
    t = tic;
    for i = 1:ntot
        % Only update display every 10k iterations
        if mod(i,10000)==0
            disp(['  ',num2str(ntot-i),' left (',...
                  datestr(toc(t)*(ntot-i)/(i*60*60*24),'DD:HH:MM:SS'),')'])
        end
        % Grab local region:
        [ii,jj,kk] = ind2sub(d(1:3),pp.ind(i));
        wimg = pp.img(max(1,ii-pp.n(1)):min(d(1),ii+pp.n(1)),...
                      max(1,jj-pp.n(2)):min(d(2),jj+pp.n(2)),...
                 	  max(1,kk-pp.n(3)):min(d(3),kk+pp.n(3)));
        wmask = pp.mask(max(1,ii-pp.n(1)):min(d(1),ii+pp.n(1)),...
                        max(1,jj-pp.n(2)):min(d(2),jj+pp.n(2)),...
                        max(1,kk-pp.n(3)):min(d(3),kk+pp.n(3)));
        % Determine percentiles:
        pthresh = pp.thresh;
        if pchk
            pthresh = prctile(wimg(wmask),pp.thresh);
        end
        for ith = 1:nth
            BW = eval(['wimg ',pp.tmode, 'pthresh(ith);']);
            BW(~wmask) = pp.defVal;
            MF(ith,1:nmf,i) = feval(func,BW,pp.voxsz);
        end
        if pchk
            MF(:,end,i) = pthresh;
        end
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



