% Calculates regional minkowski functionals using a moving window method
% Inputs:
%   img     = 2D/3D numeric/logical matrix
%   Optional Name/Value pairs:
%       thresh  = image thresholds for binarization
%       tmode   = (p indicates percentile) <,>,==,etc.
%                   defaults: > for numeric, == for logical
%       ord     = which MF to perform (default to all)
%                 ord   2D          3D
%                 0     Area        Volume
%                 1     Perimeter   Surface Area
%                 2     Euler       Mean Breadth
%                 3     n/a         Euler
%       n       = moving window radius
%       gridsp  = grid spacing
%       ind     = (3D) indices of matrix locations to analyze
%       voxsz   = voxel size in 3D
%       mask    = binary image mask
%       defVal  = (logical) value assigned to voxels outside the mask
%       prog    = logical check to display waitbar/progress
% Output:
%   MF     = [#thresh x nMF x #pts] Minkowski Functionals 
%               * 2D : Area, Perimeter, Euler-Poincare
%               * 3D : Volume, Surface Area, Mean Breadth, Euler-Poincare
%                   ** one extra MF when percentiles used (for threshold)
%   labels = cell array of MF names

function [MF,p] = minkowskiFun(img,varargin)

% Validate all inputs:
p = parseInputs(img,varargin{:});

nmf = length(p.ord);
nth = length(p.thresh);
voxvol = prod(p.voxsz);
nCth = length(p.Cthresh);
if p.calcConfl
    dCth = diff(p.Cthresh);
end

logchk = islogical(img);

% Now calculate the MF values:
if p.gchk % Perform global analysis
    
    MF = nan(nth,nmf+p.calcConfl);
    
    % Convert percentiles to threshold values:
    if p.pchk
        p.thresh = prctile(p.img(p.mask),p.thresh);
    end
    
    if p.prog, hw = waitbar(0,'Processing Minkowski Functional Analysis ...'); end
    
    % Loop over thresholds
    for ith = 1:nth
        if logchk
            BW = img;
        else
            BW = feval(p.tfun,img,p.thresh(ith));
        end
        BW(~p.mask) = p.defVal;
        if p.calcConfl
            % Smooth binary matrix for confluence calcuation:
            fprintf('Smothing binary image: win = %u ; str = %f\n',p.Cwin,p.Cstr);
            BWdens = smooth3(double(BW),'gaussian',p.Cwin,p.Cstr);
        end
        if nmf
            MF(ith,1:nmf) = feval(p.func,BW,p.voxsz,p.ord);
        end
        if p.calcConfl
            X = nan(nCth,1);
            for j = 1:nCth
                BW = feval(@gt,BWdens,p.Cthresh(j));
                X(j) = feval(p.Cfunc,BW);
            end
            X = (X-1).^2;
            MF(ith,end) = sqrt( dCth * ((X(1:end-1) + X(2:end))/2) );
        end
        if p.prog, waitbar(ith/nth,hw); end
    end
    % Adjust for voxel size (except for Euler):
    MF(:,(p.ord~=p.nD)) = MF(:,(p.ord~=p.nD));
    if p.prog, delete(hw); end
    
else % Perform analysis on moving window
    
    ntot = length(p.ind);
    disp(['Number of iterations: ',num2str(ntot)])

    pthresh = p.thresh;
    if p.pchk
        pthresh = prctile(img(p.mask),pp.thresh);
    end
        
    tt = tic;
    MF = zeros(nth,nmf+p.calcConfl,ntot);
    hw = waitbar(0,'Calculating local Minkowski Functionals ...');
    for ith = 1:nth
        fprintf(['Calculating local MF: img ',p.tmode,' %f\n'],pthresh(ith));
        if logchk
            BW = img;
        else
            BW = feval(p.tfun,img,pthresh(ith));
        end
        BW(~p.mask) = p.defVal;
        if p.calcConfl
            % Smooth binary matrix for confluence calcuation:
            fprintf('Smothing binary image: win = %u ; str = %f\n',p.Cwin,p.Cstr);
            BWdens = smooth3(double(BW),'gaussian',p.Cwin,p.Cstr);
        end
        
        % Now loop over window locations:
        for i = 1:ntot
            if i==1
                t = tic;
            end
            % Only update display every 10k iterations
            if mod(i,10000)==0
                disp(['  ',num2str(ntot-i),' left (',...
                      datestr(toc(t)*(ntot-i)/(i*60*60*24),'DD:HH:MM:SS'),')'])
            end
            
            % Grab local region:
            [ii,jj,kk] = ind2sub(p.D,p.ind(i));
            ii = max(1,ii-p.n(1)):min(p.D(1),ii+p.n(1));
            jj = max(1,jj-p.n(2)):min(p.D(2),jj+p.n(2));
            kk = max(1,kk-p.n(3)):min(p.D(3),kk+p.n(3));
            wBW = BW(ii,jj,kk);
            wmask = p.mask(ii,jj,kk);
            wvol = nnz(wmask)*voxvol*ones(1,nmf);
            wvol(p.ord==p.nD) = 1; % Don't normalize Euler value
            
            % Calculate selected MF and normalize by wvol for density:
            if nmf
                MF(ith,1:nmf,i) = feval(p.func,wBW,p.voxsz,p.ord)./wvol;
            end
            
            % Calculate confluence if desired:
            if p.calcConfl
                wBWdens = BWdens(ii,jj,kk);
                X = nan(nCth,1);
                for j = 1:nCth
                    wBW = feval(@gt,wBWdens,p.Cthresh(j));
                    X(j) = feval(p.Cfunc,wBW);
                end
                X = (X-1).^2;
                MF(ith,end,i) = sqrt( dCth * ((X(1:end-1) + X(2:end))/2) );
            end
            waitbar(i/ntot,hw,['Complete: ',num2str(i),'/',num2str(ntot)]);
        end
    end
    delete(hw);
    toc(tt);
end

% Raw 2D MF calculation
function p = calcMF2D(BW,voxsz,ord)
p = nan(1,length(ord));
if ismember(0,ord)
    p(ord==0) = imAreaEstimate(BW,voxsz);
end
if ismember(1,ord)
    p(ord==1) = imPerimeterEstimate(BW,voxsz);
end
if ismember(2,ord)
    p(ord==2) = imEuler2dEstimate(BW);
end
  
% Raw 3D MF calculation
function p = calcMF3D(BW,voxsz,ord)
p = nan(1,length(ord));
if ismember(0,ord)
    p(ord==0) = imVolumeEstimate(BW,voxsz);
end
if ismember(1,ord)
    p(ord==1) = imSurfaceEstimate(BW,voxsz);
end
if ismember(2,ord)
    p(ord==2) = imMeanBreadth(BW,voxsz);
end
if ismember(3,ord)
    p(ord==3) = imEuler3dEstimate(BW);
end

function p = parseInputs(img,varargin)
% Validate optional inputs:
if ~(isnumeric(img)||islogical(img))
    error('Invalid input image.')
end
p = inputParser;
addParameter(p,'thresh',nan,@isvector);
addParameter(p,'tmode','>',@ischar);
addParameter(p,'ord',nan,@isnumeric);
addParameter(p,'n',inf(1,3),@isvector);
addParameter(p,'gridsp',nan,@isvector);
addParameter(p,'ind',[],@isvector);
addParameter(p,'voxsz',ones(1,3),@isvector);
addParameter(p,'mask',[],@islogical);
addParameter(p,'defVal',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'prog',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'calcConfl',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'Cwin',9,@(x)isnumeric(x)&&(x>0));
addParameter(p,'Cstr',2,@(x)isnumeric(x)&&(x>0));
addParameter(p,'Cthresh',0:0.05:1,@isvector);
parse(p,varargin{:});
p = p.Results;
D = size(img);
if (length(D)>3) || (length(D)<2)
    error('Invalid image size. Must be 2D or 3D.');
end
if islogical(img)
    p.thresh = true;
    p.tmode = '==';
elseif isnan(p.thresh)
    % If no thresholds input, use 50th percentile
    p.thresh = 0.5;
    p.tmode = 'p>';
end
if any((p.ind<1) | (p.ind>prod(D)))
    warning('Invalid matrix index ... proceeding with full image.');
    p.ind = [];
end
% Global vs Local analysis:
p.gchk = any(isinf(p.n));
if p.gchk
    p.nD = length(D);
    p.ind = [];
else
    if isempty(p.ind)
        if isnan(p.gridsp)
            p.gridsp = p.n+1;
        end
        % need to determine locations:
        p.ind = getGridInd(D,p.gridsp,p.mask);
    end
    p.nD = nnz(p.n);
end
if isnan(p.ord)
    % Default to all MF
    p.ord = 0:p.nD;
elseif any(~ismember(p.ord,0:p.nD))
    error('Invalid MF value requested.');
end
if strncmp(p.tmode,'p',1)
    p.pchk = true;
    p.tmode(1) = [];
else
    p.pchk = false;
end
switch p.tmode
    case '<'
        p.tfun = @lt;
    case '<='
        p.tfun = @le;
    case '>'
        p.tfun = @gt;
    case '>='
        p.tfun = @ge;
    case '=='
        p.tfun = @eq;
    case '~='
        p.tfun = @ne;
    otherwise
        error(['Threshold function (',p.tmode,') is not valid.']);
end
if isempty(p.mask)
    p.mask = true(D);
elseif ~all(size(p.mask)==D)
    error('Mask dimensions do not fit image');
end
% Determine which MF functions to perform:
if p.nD==2
    p.func = @calcMF2D;
    p.labels = {'Area','Perimeter','Euler'};
elseif p.nD==3
    p.func = @calcMF3D;
    p.labels = {'Volume','SurfaceArea','MeanBreadth','Euler'};
else
    error('Invalid number of dimensions.');
end
p.labels = p.labels(p.ord+1);
if p.calcConfl
    p.labels{end+1} = 'alpha';
    if p.nD==2
        p.Cfunc = @imEuler2dEstimate;
    elseif p.nD==3
        p.Cfunc = @imEuler3dEstimate;
    end
end
p.D = D;

