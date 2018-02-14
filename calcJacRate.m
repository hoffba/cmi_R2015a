function dJ = calcJacRate(dX,dt,voxsz)
% Inputs:
%   dX = Deformation map: [nx,ny,nz,3]
%   dt = time interval (scalar)
%   voxsz = voxel dimensions
% Output:
%   dJ = time-derivative determinant of spatial Jacobian
%
% Theory: 
%   (from https://en.wikiversity.org/wiki/Continuum_mechanics/Time_derivatives_and_rates)
%   J = det(F)
%   Phi = deformation field
%   Fdot = d/dX(dPhi/dt)
%   d/dt(J) = J * tr( Fdot * inv(F) )

dJ = [];
[d(1),d(2),d(3),nv] = size(dX);
if (nv==3) && isscalar(dt) && isvector(voxsz)
    
    dJ = zeros(d);
    fprintf('Calculating spatial Jacobian ...\n');
    F = def2jac(dX,voxsz);
    fprintf('Calculating spatial rate Jacobian ...\n');
    FF = def2jac(dX/dt,voxsz);
    
    np = prod(d);
    n = np*9;
    f = zeros(3); % J
    ff = zeros(3);% dJ/dt
    gi = mod(1:np,round(np/20))==0;
    fprintf('Calculating Jacobian rate determinant ...\n');
    for i = 1:np
        f(:) = F(i:np:n);
        ff(:) = FF(i:np:n);
        dJ(i) = det(f) * trace( ff / f )/3;
        if gi(i)
            fprintf('%u%% Complete\n',sum(gi(1:i))*5);
        end
    end
end