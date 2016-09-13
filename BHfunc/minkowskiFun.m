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
%       n       = [1 x nD] moving window radius
%       gridsp  = grid spacing
%       ind     = (3D) indices of matrix locations to analyze
%       voxsz   = voxel size in 3D
%       mask    = binary image mask
%       defVal  = (logical) value assigned to voxels outside the mask
%       prog    = logical check to display waitbar/progress
%       calcConfl = logical check to calculate confluence
%       Cwin    = confluence smoothing window (default 9)
%       Cstr    = confluence smoothing strength (default 2)
%       Cthresh = confluence thresholds to integrate (default 0:0.05:1)
%       gchk    = TF flag for processing global measurements (default TRUE)
% Output:
%   p = structure containing all parameters and results
%    .gMF = [#thresh x nMF] global MF
%    .MF  = [#thresh x nMF x (1+#pts)] Minkowski Functionals
%    .labels = names of MF requested
% Minkowski Functionals:
%    * First analysis is always global
%    * 2D : Area, Perimeter, Euler-Poincare
%    * 3D : Volume, Surface Area, Mean Breadth, Euler-Poincare
%      ** one extra MF when percentiles used (for threshold)

function p = minkowskiFun(img,varargin)

% Validate all inputs:
p = parseInputs(img,varargin{:});

nmf = length(p.ord);
nth = length(p.thresh);
ni = length(p.ind) + p.gchk;
voxvol = prod(p.voxsz);
nCth = length(p.Cthresh);
if p.calcConfl
    dCth = diff(p.Cthresh);
end

logchk = islogical(img);

% Initialize results matrix
p.MF = nan(nth,nmf+p.calcConfl,ni);

% Convert percentiles to threshold values:
if p.pchk
    p.thresh = prctile(p.img(p.mask),p.thresh);
end

if p.prog, hw = waitbar(0,'Processing Minkowski Functional Analysis ...'); end
fprintf('Number of iterations: %u\n',ni)

% Loop over thresholds
p.wnum = nan(ni,1); % for saving window volumes (used in normalization)
for ith = 1:nth
    fprintf('Processing threshold: %s %f\n',p.tmode,p.thresh(ith));
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
    
    % Loop over desired windows:
    for i = 1:ni
        
        if i==1
            t = tic;
        end
        % Only update display every 10k iterations
        if mod(i,10000)==0
            disp(['  ',num2str(ni-i),' left (',...
                  datestr(toc(t)*(ni-i)/(i*86400),'DD:HH:MM:SS'),')'])
        end
            
        % Determine matrices to evaluate:
        % * If global analysis is requested, it's always done first and
        %   separated out after the windowed analysis.
        if (i==1) && p.gchk % GLOBAL
            
            wBW = BW;
            wmask = p.mask;
            if p.calcConfl
                wBWdens = BWdens;
            end
            
        else % LOCAL GRIDDED WINDOW
            
            % Grab local region:
            [ii,jj,kk] = ind2sub(p.D,p.ind(i-p.gchk));
            ii = max(1,ii-p.n(1)):min(p.D(1),ii+p.n(1));
            jj = max(1,jj-p.n(2)):min(p.D(2),jj+p.n(2));
            kk = max(1,kk-p.n(3)):min(p.D(3),kk+p.n(3));
            wBW = BW(ii,jj,kk);
            wmask = p.mask(ii,jj,kk);
            if p.calcConfl
                wBWdens = BWdens(ii,jj,kk);
            end
            
        end
        
        % Calculate MF normalization values:
        if isnan(p.wnum(i))
            p.wnum(i) = nnz(wmask);
        end
        normval = p.wnum(i)*voxvol*ones(1,nmf); % For normalizing V, S, B
        normval(p.ord==p.nD) = p.wnum(i); % For normalizing X and alpha
        
        % Calculate MF values:
        if nmf
            p.MF(ith,1:nmf,i) = feval(p.func,wBW,p.voxsz,p.ord)./normval;
        end
        
        % Calculate confluence:
        if p.calcConfl
            X = nan(nCth,1);
            for j = 1:nCth
                cBW = feval(@gt,wBWdens,p.Cthresh(j));
                X(j) = feval(p.Cfunc,cBW);
            end
            X = (X-1).^2;
            p.MF(ith,end,i) = sqrt( dCth * ((X(1:end-1) + X(2:end))/2) )/p.wnum(i);
        end
        
    end
    
    if p.prog, waitbar(ith/nth,hw); end
end

% Separate out the GLOBAL results:
if p.gchk
    p.gMF = p.MF(:,:,1);
    p.gnum = p.wnum(1);
    p.wnum(1) = [];
    p.MF(:,:,1) = [];
end

if p.prog, delete(hw); end

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
% Validate inputs:
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
addParameter(p,'gchk',true,@(x)islogical(x)&&isscalar(x));
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
if any(isinf(p.n))
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

