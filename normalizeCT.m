function Inorm = normalizeCT(I,BW,voxsz)
% Normalization function based on Gallardo-Estrella et al. 2016

n = 6;
s = 2.^(0:(n-2));
r = [70.34, 67.54, 60.90, 51.45, 36.14];

mI = I;
mI(~BW) = nan;

% Initialize to F^n:
% Inorm = imgaussfilt(I,s(1)./voxsz(1:2));
Inorm = filtImg(mI,genFilt(s(end)./voxsz));

% L_0 is original image:
Lprev = mI;

for i = 1:n-1 % Loop over kernels
%     L = imgaussfilt(I,s(i)./voxsz(1:2));
    L = filtImg(mI,genFilt(s(i)./voxsz));
    F = Lprev - L;
    Inorm = Inorm + r(i)/std(F(BW))*F;
    if i<n-1
        Lprev = L;
    end
end

Inorm(~BW) = I(~BW);

function h = genFilt(sig)
r = ceil(2*sig);
[X,Y] = ndgrid(-r(1):r(1),-r(2):r(2));
h = exp( -X.^2/(2*sig(1)^2)  -Y.^2/(2*sig(2)^2) ) / sqrt(((2*pi)^3 * prod(sig.^2)));
h = h/sum(h(:));

function L = filtImg(I,h)
d = size(I);
L = zeros(d);
for i = 1:d(3)
    L(:,:,i) = nanconv(I(:,:,i),h,'edge','nanout');
end