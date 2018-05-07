function [E,V] = eig3d(J) 
% Calculates the eigenvalue and vector maps of 3D full spatial Jacobian matrix
% (From Amelon et al. 2011)
% Inputs:
%   J = [d1 x d2 x d3 x 9] full Jacobian matrix
%           (4D in order of xx, xy, xz, yx, yy, yz, zx, zy, zz)
% Outputs:
%   D = [d1 x d2 x d3 x 3] maps of sorted eigenvalues
%   V = [d1 x d2 x d3 x 3 x 3] maps of sorted corresponding eigenvectors

E = [];
V = [];
[d(1),d(2),d(3),nJ] = size(J);

if nJ==9
    
    % Factor out rotational tensor (makes it symmetrical):
    % (Amelon et al. 2011)
    % F = J'*J
    % F = [ xx, xy, xz, yy, yz, zz ]
    disp('Factoring out rotational tensor ...')
    F(:,:,:,1) = sum(J(:,:,:,[1,4,7]).^2,4);
    F(:,:,:,2) = sum(J(:,:,:,[1,4,7]).*J(:,:,:,[2,5,8]),4);
    F(:,:,:,3) = sum(J(:,:,:,[1,4,7]).*J(:,:,:,[3,6,9]),4);
    F(:,:,:,4) = sum(J(:,:,:,[2,5,8]).^2,4);
    F(:,:,:,5) = sum(J(:,:,:,[2,5,8]).*J(:,:,:,[3,6,9]),4);
    F(:,:,:,6) = sum(J(:,:,:,[3,6,9]).^2,4);
    clear J
        
    % Derive invariants:
    disp('Deriving invariants ...')
    I1 = (F(:,:,:,1) + F(:,:,:,4) + F(:,:,:,6))/3;
    I2 = (F(:,:,:,1).*F(:,:,:,4)+F(:,:,:,1).*F(:,:,:,6)+F(:,:,:,4).*F(:,:,:,6)) ...
            -(F(:,:,:,2).^2+F(:,:,:,3).^2+F(:,:,:,5).^2);
    I3 = F(:,:,:,1).*F(:,:,:,4).*F(:,:,:,6)+2*F(:,:,:,2).*F(:,:,:,3).*F(:,:,:,5) ...
            -(F(:,:,:,6).*F(:,:,:,2).^2+F(:,:,:,4).*F(:,:,:,3).^2+F(:,:,:,1).*F(:,:,:,5).^2);
        
    % Calculate characteristics:
    disp('Calculating characteristics ...')
    s = I1.^3 - I1.*I2/2 + I3/2;
    clear I3
    v = I1.^2 - I2/3;
    clear I2
    phi = acos(s./v.*sqrt(1./v))/3;
    
    % Generate eigenvalues:
    disp('Generating eigenvalues ...')
    v = 2*sqrt(v);
    E = nan([d,3]);
    E(:,:,:,1) = real(sqrt(I1 + v.*cos(phi)));
    E(:,:,:,2) = real(sqrt(I1 - v.*cos(pi/3+phi)));
    E(:,:,:,3) = real(sqrt(I1 - v.*cos(pi/3-phi)));
    clear I1
    
    % Generate eigenvectors:
    if nargout==2        
        disp('Generating eigenvectors ...')
        V = nan([d,3,3]);
        for i = 1:2
            A = F(:,:,:,1) - E(:,:,:,i);
            B = F(:,:,:,4) - E(:,:,:,i);
            C = F(:,:,:,6) - E(:,:,:,i);
            V(:,:,:,i,1) = (F(:,:,:,2).*F(:,:,:,5) - B.*F(:,:,:,3)) ...
                .* (F(:,:,:,3).*F(:,:,:,5) - C.*F(:,:,:,1));
            V(:,:,:,i,2) = (F(:,:,:,3).*F(:,:,:,5) - C.*F(:,:,:,2)) ...
                .* (F(:,:,:,3).*F(:,:,:,2) - A.*F(:,:,:,5));
            V(:,:,:,i,3) = (F(:,:,:,2).*F(:,:,:,5) - B.*F(:,:,:,3)) ...
                .* (F(:,:,:,3).*F(:,:,:,2) - A.*F(:,:,:,5));
        end
        V(:,:,:,3,1) = V(:,:,:,1,2).*V(:,:,:,2,3) - V(:,:,:,2,2).*V(:,:,:,1,3);
        V(:,:,:,3,2) = V(:,:,:,1,3).*V(:,:,:,2,1) - V(:,:,:,1,1).*V(:,:,:,2,3);
        V(:,:,:,3,3) = V(:,:,:,1,1).*V(:,:,:,2,2) - V(:,:,:,1,2).*V(:,:,:,2,1);
    end
    
    disp('... done')
end