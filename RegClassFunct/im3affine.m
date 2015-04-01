function imout = im3affine(imin, T, method)
%%
% Inputs:
%   imin          - 3D image volume
%   T             - 4x4 affine transform matrix
%   method        - interpolation method (1=trilinear,2=nearest neighbor)
%%
if ndims(imin) == 3
    if ~isequal(size(T),[4 4])
        errordlg('Transform should be a 4x4 affine matrix for 3D volume.');
    end
else
    errordlg('Input image should be 3D volume.');
end
if ~exist('method','var') || isempty(method) || any(strcmp(method,{'linear','nearest'}))
   method = 'linear';
end

d=size(imin);
n = numel(imin);
[X,Y,Z] = meshgrid(0:d(1)-1,0:d(2)-1,0:d(3)-1);
C=[reshape(X,[1,n]);reshape(Y,[1,n]);reshape(Z,[1,n]);ones(1,n)];
Ci=inv(T) * C;
Xi=reshape(Ci(1,:),d);
Yi=reshape(Ci(2,:),d);
Zi=reshape(Ci(3,:),d);
imout = interp3(X,Y,Z,imin,Xi,Yi,Zi);




