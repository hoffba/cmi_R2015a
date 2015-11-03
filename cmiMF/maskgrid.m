
function ind = maskgrid(mask,r)
% Find grid indices on mask for sub-sampling
% Inputs:
%   mask = binary image mask
%   r = vector of grid spacing integers

if nargin==2 && islogical(mask) && (length(r)==3)
    
    % Find limits of VOI:
    mchk = isempty(mask);
    if mchk
        mmin = [ find(max(max(mask,[],2),[],3),1) ,...
                 find(max(max(mask,[],1),[],3),1) ,...
                 find(max(max(mask,[],1),[],2),1) ];
    else
        mmin = ones(1,3);
    end
    
    % Generate grid:
    gridmask = false(self.dims(1:3));
    gridmask( mmin(1):pp.n(1):end , ...
              mmin(2):pp.n(2):end , ...
              mmin(3):pp.n(3):end ) = true;
    if mchk
        gridmask = gridmask & mask;
    end
    ind = find(gridmask);
    
end

