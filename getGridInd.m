% Determine grid indices within matrix
% Inputs:
%       d    = matrix extents
%       gsp  = grid spacing [1x3] vector or scalar for iso-grid
%       mask = (optional) binary mask
% Output:
%       ind = vector of matrix indices
function ind = getGridInd(d,gsp,mask)

if nargin>1
    
    if length(gsp)==1
        gsp = gsp*ones(1,3);
    end
    
    % Determine min extent of mask:
    if (nargin<3)||isempty(mask)
        mmin = ones(1,3);
        dm = nan;
    else
        dm = size(mask);
        if ~all(d==dm)
            error('Mask size does not equal input dimensions.');
        else
            mmin = [ find(max(max(mask,[],2),[],3),1) ,...
                     find(max(max(mask,[],1),[],3),1) ,...
                     find(max(max(mask,[],1),[],2),1) ];
        end
    end
    
    [X,Y,Z] = meshgrid(mmin(2):gsp(2):d(2),...
                       mmin(1):gsp(1):d(1),...
                       mmin(3):gsp(3):d(3));
    ind = sub2ind(d,Y(:),X(:),Z(:));
    if ~isnan(dm)
        ind(~mask(ind)) = [];
    end
    
end