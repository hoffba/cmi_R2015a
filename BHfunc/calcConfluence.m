function c = calcConfluence(I,varargin)
% c = calcConfluence(I)
% Calculates confluence of image I
%   - if I is logical, blur first for density
% Inputs:
%   I = image matrix
%   Name/Value pairs:
%       thresh  = thresholds to assess
%       n       = moving window radius
%       ind     = (3D) indices of matrix locations to analyze
%       mask    = binary image mask
%       defVal  = (logical) value assigned to voxels outside the mask
%       prog    = logical check to display waitbar/progress
% Output:
%   c = confluence value [1 x np]

p = inputParser;
addRequired(p,'I',@(x)isnumeric(x)||islogical(x));
addParameter(p,'thresh',[],@isvector);
addParameter(p,'tmode','>',@ischar);
addParameter(p,'n',inf(1,3),@isvector);
addParameter(p,'ind',[],@isvector);
addParameter(p,'mask',[],@islogical);
addParameter(p,'defVal',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'prog',false,@(x)islogical(x)&&isscalar(x));
parse(p,I,varargin{:});
thresh = p.Results.thresh;
tmode = p.Results.tmode;
n = p.Results.n;
ind = p.Results.ind;
mask = p.Results.mask;
defVal = p.Results.defVal;
prog = p.Results.prog;
clear p;

if isempty(mask)
    mask = true(size(I));
elseif any(size(mask)~=size(I))
    error('Mask and Image dimensions must match');
end
switch tmode
    case '<'
        tfun = @lt;
    case '<='
        tfun = @le;
    case '>'
        tfun = @gt;
    case '>='
        tfun = @ge;
    case '=='
        tfun = @eq;
    case '~='
        tfun = @ne;
    otherwise
        error(['Threshold function (',tmode,') is not valid.']);
end
if ~ismember(tmode,{'<','<=','>','>=','==','~='})
    error(['Threshold function (',tmode,') is not valid.']);
end

% Blur the binary image to get density map:
if islogical(I)
    str = 2; % 3D Gaussian StDev in voxels
    win = min(min(size(I)),round(4*str+1));
    fprintf('Smothing binary image: win = %u ; str = %f\n',win,str);
    I = smooth3(double(I),'gaussian',win,str);
end

% Automatically determine thresholds
if isempty(thresh)
    Imin = min(I(:));
    Imax = max(I(:));
    thresh = Imin:(Imax-Imin)/10:Imax;
elseif length(thresh)==1
    error('Number of thresholds must be >1.');
end

if isempty(n) || any(isinf(n))
    gchk = true;
else gchk = false;
end

if ~gchk && isempty(ind)
    ind = find(mask);
end

nth = length(thresh);
np = max(length(ind),1);
X = zeros(nth,np);
if gchk % Global
    
    % Calculate Euler at multiple density thresholds
    if prog, hw = waitbar(0,'Calculating Euler over density thresholds ...'); end
    for i = 1:nth
        BW = feval(tfun,I,thresh(i));
        BW(~mask) = defVal;
        X(i) = imEuler3dEstimate(BW);
        if prog, waitbar(i/nth,hw); end
    end
    if prog, delete(hw); end

else % Moving window
    
    d = size(I);
    disp(['Number of iterations: ',num2str(np)])

    if prog, hw = waitbar(0,'Calculating local Minkowski Functionals ...'); else hw = []; end
    t = tic;
    for i = 1:np
        % Only update display every 10k iterations
        if mod(i,10000)==0
            disp(['  ',num2str(np-i),' left (',...
                  datestr(toc(t)*(np-i)/(i*60*60*24),'DD:HH:MM:SS'),')'])
        end
        % Grab local region:
        [ii,jj,kk] = ind2sub(d(1:3),ind(i));
        ii = max(1,ii-n(1)):min(d(1),ii+n(1));
        jj = max(1,jj-n(2)):min(d(2),jj+n(2));
        kk = max(1,kk-n(3)):min(d(3),kk+n(3));
        wimg = I(ii,jj,kk);
        wmask = mask(ii,jj,kk);
        nm = nnz(wmask);
        for ith = 1:nth
            BW = feval(tfun,wimg,thresh(ith));
            BW(~wmask) = defVal;
            X(ith,i) = imEuler3dEstimate(BW)/nm;
        end
        if prog, waitbar(i/np,hw,['Complete: ',num2str(i),'/',num2str(np)]); end
    end
    if prog, delete(hw); end
    toc(t);
    
end

% Calculate confluence (clumpiness) using trapezoidal integration:
X = X.^2;
c = sqrt( diff(thresh(:))' * ((X(1:end-1,:) + X(2:end,:))/2) );

