function [T1,a,b,res] = T1map_mIR(img,ti,mask)
% Model: S = a + b*exp(-TI/T1)
% Finds estimates of T1, a, and b using nonlinear least squares
% Inputs: 
%   img = [ny nx nz nTI] image matrix
%   ti  = [1 nTI] inversion times
% Outputs:
%   T1  = 3D T1 map
%   a   = 3D a map
%   b   = 3D b map
%   res = 3D map of rms error for the fit

[d(1),d(2),d(3),ni] = size(img);

if nargin<3
    mask = true(1);
end

if isreal(img)
    func = 'rdNlsPr';
else
    func = 'rdNls';
end

