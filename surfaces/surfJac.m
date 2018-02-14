function h = surfJac(N,mask,voxsz,jac,clim,dt)
% Generate surface from mask with Jacobian map overlaid

nd = size(jac,4);
if nargin<6
    dt = 1;
end

if nd==1
    
    % 3D Jacobian determinant:
    h = surfmap(N,mask,voxsz,jac.^(1/dt),clim);
%     savesurf(h,1:2);
    
elseif nd==9
    
    % 2D Jacobian determinant tangent to surface:
    [fv,x,y,z] = mask2surf(mask,voxsz);
    
    % Interpolate Jacobian at vertices:
    nv = length(fv.vertices);
    ijac = zeros(nv,9);
    for i = 1:9
        ijac(:,i) = interp3(x,y,z,jac(:,:,:,i),fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3));
    end
    ijac = permute(reshape(ijac',3,3,[]),[2,1,3]);

    % Determine planar Jacobian from isonormals:
    n = isonormals(x,y,z,mask,fv.vertices);
    pjac = zeros(nv,1);
    for i = 1:size(n,1)
        try
        R = vrrotvec2mat(vrrotvec([0,1,0],n(i,:)));
        catch err
            disp(num2str(i))
        end
        tjac = R*ijac(:,:,i)*R';
        pjac(i) = det(tjac(1:2,1:2));
    end
    pjac = pjac.^(1/dt);
    
    % Generate surface figure:
    h = surfmap(N,fv,pjac,clim);
    
    % Add normal vectors:
%     isub = round(linspace(1,length(n),1000));
%     hold(h.axes1,'on');
%     quiver3(h.axes1,fv.vertices(isub,1),fv.vertices(isub,2),fv.vertices(isub,3),...
%         n(isub,1),n(isub,2),n(isub,3));
%     hold(h.axes1,'off');
    
else
    error('Invalid Jacobian input.');
end




