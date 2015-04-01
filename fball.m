function H = fball(r)
% Creates 3D ball filter (or circle if 2D)
if (nargin==2) && isnumeric(r)
    nd = length(r);
    if nd>1
        r = round(r);
        dim = -r:r;
        if nd==2
            [X,Y] = ndgrid(dim,dim);
            Z = 0;
        elseif nd==3
            [X,Y,Z] = ndgrid(dim,dim,dim);
        end
        H = double((X.^2 + Y.^2 + Z.^2) <= (1.05*r)^2);
    end
end
    