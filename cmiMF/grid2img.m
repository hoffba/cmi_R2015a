
function img = grid2img(MF,ind,mask)
% Interpolate grid data to original dimensions
% Inputs:
%   MF = [nx1] vector of values
%   ind = [nx1] vector matrix indices
%   mask = binary mask 3D matrix

if nargin==3 && (length(MF)==length(ind)) && all((ind>0)&(ind<=numel(mask)))
    
    % Generate interpolation model:
    d = size(mask);
    [Xq,Yq,Zq] = meshgrid(1:d(2),1:d(1),1:d(3));
    [Y,X,Z] = ind2sub(d,ind);
    F = scatteredInterpolant(X,Y,Z,MF(:),'linear','none');
    
    % Calculate full grid
    img = F(Xq,Yq,Zq);
    
    % Validate results:
    img(~mask | isnan(img)) = 0;
    
end



