% Affine spatial coregistration
% affineReg(Im0,Vox0,Im1,Vox1,metric,opts)
% Inputs:   Im0 = Static image
%           Vox0 = Voxel dimensions of Im0
%           Im1 = Moving image
%           Vox1 = Voxel dimensions of Im1
%           affType = Definition of affine terms used for transform
%               'RT' = Rotate / Translate
%               'RTiS' = Rotate / Translate / isotropic Scale
%               'RTS' = Rotate / Translate / Scale
%               'RTSK' = Rotate / Translate / Scale / Shear
%           metric = Cost function of optimization
%               'mi' = Mutual Information
%               'sd' = Squared Difference
%           interpm = Interpolation method
%               0: linear interpolation and outside pixels set to nearest pixel
%               1: linear interpolation and outside pixels set to zero
%               2: cubic interpolation and outsite pixels set to nearest pixel
%               3: cubic interpolation and outside pixels set to zero
%           opts = Structure of optimization options or optimset
function [ImReg,M] = affineReg(Im0,Im1,Vox0,Vox1,Mask0,Mask1,affType,metric,interpm,opts)

if (nargin>1) && ~(isempty(Im0) || isempty(Im1))
    d0 = size(Im0);
    d1 = size(Im1);
    chk3D = (length(d0)==3) || (length(d1)==3);
    
    % Set defaults for non-input variables
    if ~exist('Vox0','var') || isempty(Vox0)
        Vox0 = ones(1,nd);
    end
    if ~exist('Vox1','var') || isempty(Vox1)
        Vox1 = ones(1,nd);
    end
    if ~exist('Mask0','var') || isempty(Mask0)
        Mask0 = true;
    end
    if ~exist('Mask1','var') || isempty(Mask1)
        Mask1 = true;
    end
    types = {'RT','RTiS','RTS','RTSK'};
    if ~exist('affType','var') || isempty(affType)
        affType = 'RTS';
    elseif ~any(strcmp(affType,types))
        error('Invalid affine transform type input, options = RT, RTiS, RTS, or RTSK');
    end
    costfuncs = {'mi','sd'};
    if ~exist('metric','var') || isempty(metric) || ~ischar(metric)
        metric = 'mi';
    elseif ~any(strcmp(metric,costfuncs))
        error('Invalid cost function input, options = mi or sd');
    end
    if ~exist('opts','var') || isempty(opts)
        opts = struct('GradObj','off',...
                      'Display','iter',...
                      'MaxIter',100,...
                      'HessUpdate','lbfgs',...
                      'TolFun',10^-3,...
                      'TolX',10^-3);
    end

    % Initialize affine matrix to identity:
    switch affType
        case 'RT'
            if chk3D
                %    [Tx,Ty,Tz,Rx,Ry,Rz]
                x0 = [0 0 0 0 0 0];
            else
                %    [Tx,Ty,R]
                x0 = [0 0 0];
            end
        case 'RTiS'
            if chk3D
                %    [Tx,Ty,Tz,Rx,Ry,Rz,S]
                x0 = [0 0 0 0 0 0 1];
            else
                %    [Tx,Ty,R,S]
                x0 = [0 0 0 1];
            end
        case 'RTS'
            if chk3D
                %    [Tx,Ty,Tz,Rx,Ry,Rz,Sx,Sy,Sz]
                x0 = [0 0 0 0 0 0 1 1 1];
            else
                %    [Tx,Ty,R,Sx,Sy]
                x0 = [0 0 0 1 1];
            end
        case 'RTSK'
            if chk3D
                %    [Tx,Ty,Tz,Rx,Ry,Rz,Sx,Sy,Sz,Kxy,Kxz,Kyx,Kyz,Kzx,Kzy]
                x0 = [0 0 0 0 0 0 1 1 1 0 0 0 0 0 0];
            else
                %    [Tx,Ty,R,Sx,Sy,Kxy,Kyx]
                x0 = [0 0 0 1 1 0 0];
            end
    end
    
    % Determine image space (origin at center of 3D image):
        v0 = Vox0 .* (1 - d0)/2;
        v1 = Vox0 .* (d0 - 1)/2;
        [xx,yy,zz] = meshgrid(linspace(v0(2),v1(2),d0(2)),...
                              linspace(v0(1),v1(1),d0(1)),...
                              linspace(v0(3),v1(3),d0(3)));
    xyz0 = cat(4,xx,yy,zz);
        v0 = Vox1 .* (1 - d1)/2;
        v1 = Vox1 .* (d1 - 1)/2;
        [xx,yy,zz] = meshgrid(linspace(v0(2),v1(2),d1(2)),...
                              linspace(v0(1),v1(1),d1(1)),...
                              linspace(v0(3),v1(3),d1(3)));
    xyz1 = [xx(:),yy(:),zz(:)]';
        clear xx yy zz v0 v1
    
    % Perform optimization:
    x = fminlbfgs(@(x)fitFun(x,Im0,Im1,d0,d1,xyz0,xyz1,Mask0,Mask1,metric,chk3D),...
                  x0,opts);
    
    % Perform final transformation using optimized matrix:
    M = getAffTform(x,chk3D);
    xyzReg = reshape(M * xyz1,[d1,3]);
    ImReg = interpn(xyzReg(:,:,:,1),xyzReg(:,:,:,2),xyzReg(:,:,:,3),Im1,...
                    xyz0(:,:,:,1),xyz0(:,:,:,2),xyz0(:,:,:,3));
end

% Fitting function for optimizer
function e = fitFun(x,Im0,Im1,d0,d1,xyz0,xyz1,Mask0,Mask1,metric,chk3D)
% Determine affine inputs:
nx = length(x);
s = [];
k = [];
if chk3D
    t = x(1:3);
    r = x(4:6);
    if nx>6
        s = x(7:min(nx,9));
    end
    if nx==15
        k = x(10:15);
    end
else
    t = x(1:2);
    r = x(3);
    if nx>3l
        s = x(4:min(nx,5));
    end
    if nx==7
        k = x(6:7);
    end
end
M = makeAffTform(r,t,s,k);
% ImReg = affine_transform_3d_double(Im1,M,mode);
xyzReg = reshape(M * xyz1,[d1,3]);
ImReg = interpn(xyzReg(:,:,:,1),xyzReg(:,:,:,2),xyzReg(:,:,:,3),Im1,...
                xyz0(:,:,:,1),xyz0(:,:,:,2),xyz0(:,:,:,3));
if numel(Mask1)>1
    % MaskReg = logical(round(affine_transform_3d_single(single(Mask1),single(M),single(0))));
    MaskReg = interpn(xyzReg(:,:,:,1),xyzReg(:,:,:,2),xyzReg(:,:,:,3),Mask1,...
                    xyz0(:,:,:,1),xyz0(:,:,:,2),xyz0(:,:,:,3)) & Mask0;
elseif numel(Mask0)==1
    MaskReg = true(d0);
end
e = calcCostFcn(Im0(MaskReg),ImReg(MaskReg),metric);




