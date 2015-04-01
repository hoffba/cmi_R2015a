function W = kernelNOVA(imgDiff, sdev, p, d)
% W = trapezoidKernel(imgDiff, sdev, p, d)
% Weighting function for NOVA filter
% Inputs: imgDiff = (I - Io) where Io is the window center value
%                   * can be single value or matrix
%         sdev    = local noise in image
%                   * can be single value or matrix same size as imgDiff
%         p       = trapezoidal plateau (positive real number)
%         d       = trapezoidal extent (positive real number)
% Output: W       = resulting weight (between 0 and 1)
%
% 1_|         __________
%   |        /          \
%   |       /            \
%   |      /              \
% 0_|_____/                \______
%       -d  -p    0   p   d

m = -1/(d-p);
b = -m*d;
W = m*abs(imgDiff./sdev) + b;
W(W>1) = 1;
W(W<0) = 0;