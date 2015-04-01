% Warp 3D image using control points
function dimg = cmi_warp(img,vi,ipts,dpts,dd,vd)
% Inputs:   img     = 2D or 3D image matrix
%           vi      = [1 x (2 or 3)] vector of input image voxel dimensions
%           ipts    = [n x (2 or 3)] matrix of image points - [y,x,(z)]
%           dpts    = [n x (2 or 3)] matrix of destination points
%           dd      = [1 x (2 or 3)] vector of destination dimensions
%           vd      = [1 x (2 or 3)] vector of destination image voxel dimensions

di = size(img); nd = length(di);
np = size(ipts,1);
if (nargin==6) && (size(ipts,2)==nd) && (size(dpts,2)==nd) ...
        && (length(dd)==nd) && (size(dpts,1)==np) && any(nd==[2,3]) ...
        && (length(vi)==nd) && (length(vd)==nd)
    if nd==2
        di(3) = 1;
        dd(3) = 1;
        vi(3) = 1;
        ipts(:,3) = ones(np,1);
        dpts(:,3) = ones(np,1);
    end
    [yi,xi,zi] = meshgrid( (0.5:di(1))*vi(1) ,...
                           (0.5:di(2))*vi(2) ,...
                           (0.5:di(3))*vi(3));
    [yd,xd,zd] = meshgrid( (0.5:dd(1))*vd(1) ,...
                           (0.5:dd(2))*vd(2) ,...
                           (0.5:dd(3))*vd(3));
    ptsR = TPS3D(dpts,ipts,[yd(:),xd(:),zd(:)]);
    yR = reshape(ptsR(:,1),dd);
    xR = reshape(ptsR(:,2),dd);
    zR = reshape(ptsR(:,3),dd);
    dimg = interp3(yi,xi,zi,img,yR,xR,zR);
end