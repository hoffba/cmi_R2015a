function [A,label,fov] = jacmat2shapes(fullname)
% Loads full Jacobian data (MHD) from file, then performs eigenvalue
% analysis and calculates shape indices.
% Input:
%       fname = full file name of MHD containing full Jacobian
% Output:
%       A = 4D map containing: 

A = [];
label = {};
fov = [];

if nargin==0
    [fname,fpath] = uigetfile('*.mhd','Select fullJacobian MHD:');
    if ischar(fullname)
        fullname = fullfile(fpath,fullname);
    else
        fullname = [];
    end
else
    [fpath,fname] = fileparts(fullname);
end

if ischar(fullname) && exist(fullname,'file') && strcmp(fullname(end-3:end),'.mhd')
    disp('Loading jacobian maps...');
    [J,~,fov] = readMHD(fullname);
    d = size(J);
    if d(4)==9
        % Calculate eigenvalues/vectors:
        E = eig3d(J);
        saveMHD(fullfile(fpath,[fname(8:end),'.mhd']),E,{'E1','E2','E3'},fov);
        % Calculate shape indices:
%         label = {'Det','CL','CP','CS','ADI','FA'};
        label = {'ADI','FA'};
        A = tensorShape(E,label);
        saveMHD(fullfile(fpath,[fname(8:end),'.mhd']),A,label,fov);
    end
end