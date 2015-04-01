function H = gauss3D(s,r)
% Creates a 3D Gaussian filter
% Inputs: s = standard deviation (voxels)
%         r = (optional) radius of filter kernel

if (nargin>0) && isnumeric(s)
    if (nargin<2) || ~isnumeric(r)
        r = 5*ceil(s);
    else
        r = round(r);
    end
    dim = (-r:r);
    [X,Y,Z] = ndgrid(dim,dim,dim);
    H = exp(-(X.^2+Y.^2+Z.^2)/(2*s^2)) / ((2*pi)^(3/2) * s^3);
else
    H = [];
end