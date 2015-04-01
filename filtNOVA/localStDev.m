function sDevMap = localStDev(img,n)
% sDevMap = localStDev(img, n)
% Calculates local standard deviation map using moving window with radius n
% Inputs: img = 3D image matrix [iy,ix,iz]
%         n   = radius of moving window [2n+1,2n+1,2n+1]
% Output: sDevMatp = map of local standard deviations, excluding edges
%                    

d = size(img)-2*n;
N = (2*n+1)^3;
yy = (1:d(1))+n;
xx = (1:d(2))+n;
zz = (1:d(3))+n;

x = zeros(d);
x2 = x;
for i = -n:n
    for j = -n:n
        for k = -n:n
            x = x + img(i+yy,j+xx,k+zz);
            x2 = x2 + img(i+yy,j+xx,k+zz).^2;
        end
    end
end
sDevMap = sqrt( (x2 - 2*(x.^2)/N + (x.^2)/N) / (N-1) );
