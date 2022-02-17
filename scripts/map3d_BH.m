function m = map3d_BH(func,img,mask,varargin)
% map3d_BH applies functions to local image moving window.
%   map3d_BH(func, img, mask, 'Name', Value,...) applies the functions in
%   'func' to the input image 'img' on a local moving window within the 
%   masked region ('mask') determined by grid spacing and window radius 
%   (cubic), both determined as voxels. Values are then interpolated back
%   to the image dimensions to comprise parameter maps.
%
%   Input:  
%       (Required):
%           -func:      Cell array of cell arrays containing function
%                       handle and additional inputs.
%                   e.g. { {@(img,mask)mean(img(mask))};...
%                          {@(img,mask)mean(img,'all')} };
%           -img:       Image to be processed
%           -mask:      Logical matrix identifying region to process
%       (Optional Name/Value pairs):
%           'winRadius':    Radius of cubic window (voxels, default [10,10,10])
%           'winStep':      Increment between window locations (voxels)
%
%   The following must be adhered to:
%
%       - 'img' and 'mask' must be 3D and agree in size
%       - Function handles must include the two inputs 'img' and 'mask'
%           followed by any additional input variables. Additional inputs
%           are listed in the cell array after the function handle.
%       - The length of windowRad must equal the dimensionality of the image
%
%   Ben Hoff - 2/17/2022

if ~iscell(func)
    warning('First input must be cell array of cells containing function handle and extra inputs.');return;
end

tStart = tic;

% Parse inputs
d = size(mask,1,2,3);
p = parseInputs(varargin); % Next parse windowing options (Name/Value pairs)

% Find grid indices
ind = getGridInd(d,p.winStep,mask);
ni = numel(ind);

% Initialize function-related parameters
Nfunc = numel(func);
Nin = zeros(Nfunc,1);
win_var = cell(Nfunc,1);
win_flag = cell(Nfunc,1);
for i = 1:Nfunc
    Nin(i) = numel(func{i})-1;
    win_var{i} = cell(1,Nin(i));
    % inputs with same size as image matrix will be windowed:
    win_flag{i} = cellfun(@(x)all(size(x,1,2,3)==d),func{i}(2:end));
    win_var{i}(~win_flag{i}) = func{i}([false,~win_flag{i}]);
end

% Initialize loop values
vals = nan(ni,Nfunc); % must be in 3rd dimension for grid2img function
dispiter = round(ni/100);

h = waitbar(0,'3Dmap Calculating ...');
for i = 1:ni
    
    % Update waitbar
    if ~mod(i,dispiter)
        waitbar(i/ni, h, sprintf('3Dmap calc: %2.0f%%',100*i/ni));
    end
    
    % Extract window
    [ii,jj,kk] = ind2sub(d,ind(i));
    ii = max(1,ii-p.winRadius(1)) : min(d(1),ii+p.winRadius(1));
    jj = max(1,jj-p.winRadius(2)) : min(d(2),jj+p.winRadius(2));
    kk = max(1,kk-p.winRadius(3)) : min(d(3),kk+p.winRadius(3));
    wimg = img(ii,jj,kk);
    wmask = mask(ii,jj,kk);
    
    % Loop over functions:
    for j = 1:Nfunc
        
        % Window additional inputs
        for k = 1:Nin(j)
            if win_flag{j}(k)
                win_var{j}(k) = func{j}{k+1}(ii,jj,kk);
            end
        end
        
        % Evaluate function
        vals(i,j) = feval(func{j}{1},wimg,wmask,win_var{j}{:});
    end
    
end
delete(h);

% Map values to original dimensions
m = grid2img(permute(vals,[3,2,1]),ind,mask,3,1);

tEnd = toc(tStart);
fprintf('Total time : %4.1f min\n',tEnd/60)



function p = parseInputs(V)

p = inputParser;
addParameter(p,'winRadius',[10,10,10],@(x)isnumeric(x)&&ismember(numel(x),[1,3]));
addParameter(p,'winStep',[5,5,5],@(x)isnumeric(x)&&ismember(numel(x),[1,3]));
parse(p,V{:});
p = p.Results;

if isscalar(p.winRadius)
    p.winRadius = p.winRadius*ones(1,3);
end
if isscalar(p.winStep)
    p.winStep = p.winStep*ones(1,3);
end

