function [r,t,s,k] = decompAffine(M)
% Decomposes affine transform matrix into rotate/translate/scale/shear

t = M(1:3,end)';
M(1:3,end) = 0;
