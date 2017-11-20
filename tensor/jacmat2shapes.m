function [A,label,voxsz] = jacmat2shapes(fname)
% Loads full Jacobian data (MHD) from file, then performs eigenvalue
% analysis and calculates shape indices.
% Input:
%       fname = full file name of MHD containing full Jacobian
% Output:
%       A = 4D map containing: 

A = [];
label = {};
voxsz = [];

if nargin==0
    [fname,fpath] = uigetfile('*.mhd','Select fullJacobian MHD:');
    if ischar(fname)
        fname = fullfile(fpath,fname);
    else
        fname = [];
    end
end

if ischar(fname) && exist(fname,'file') && strcmp(fname(end-3:end),'.mhd')
    disp('Loading jacobian maps...');
    [J,~,fov] = readMHD(fname);
    d = size(J);
    voxsz = fov./d(1:3);
    if d(4)==9
        % Calculate eigenvalues/vectors:
        E = eig3d(J);
        % Calculate shape indices:
        label = {'Det','CL','CP','CS'};
        A = tensorShape(E,label);
    end
end