% Generate 3D ellipsoid structure element
function se = bwellipsoid(r)
% Input: r = [1 x 3] vector of radii
r = r(:)';
% SE matrix size
d = 2*r+1;
% SE matrix center
if r(3)==0
    ze = 0;
else
    ze = -1:1/r(3):1;
end
% Create ellipsoid binary filter
[xe,ye,ze] = meshgrid( -1:1/r(1):1 ,...
                       -1:1/r(2):1 ,...
                       ze );
se = zeros(d);
se(sqrt(xe.^2 + ye.^2 + ze.^2) <= 1.1) = 1;