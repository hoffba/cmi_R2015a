function X = cmiErode(X,r)
% Dilates a 3D binary mask using a square kernel
% (square kernel is separable)
if nargin && ~isempty(X) && islogical(X)
    if (nargin==2) && isnumeric(r)
        r = round(r);
    else % default to 3
        r = 3;
    end
    H = false(2*r+1,1);
    d = size(X);
    for i = 1:length(d)
        if d(1)>50
            X = reshape(conv(X(:),H,'same'),d);
        end
        X = shiftdim(X,1);
        d = size(X);
    end
end

