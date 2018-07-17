function Inorm = normalizeCT(I,BW,voxsz)
% Normalization function based on Gallardo-Estrella et al. 2016

n = 6;
s = 2.^(0:(n-2));
r = [70.34, 67.54, 60.90, 51.45, 36.14];

% Initialize to F^n:
Inorm = imgaussfilt(I,s(end)./voxsz(1:2));

% L_0 is original image:
Lprev = I;

for i = 1:n-1 % Loop over kernels
    L = imgaussfilt(I,s(i)./voxsz(1:2));
    F = Lprev - L;
    Inorm = Inorm + r(i)/std(F(BW))*F;
    if i<n-1
        Lprev = L;
    end
end
