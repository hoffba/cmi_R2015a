% Generate full affine transform from pieces
% A = genAffineMatrix(Tvals,Rvals,Svals,Kvals)
function A = genAffineMatrix(Tvals,Rvals,Svals,Kvals)
% Inputs: 
%       Tvals : [nd x 1] vector of translations ([x,y,z,...])
%       Rvals : [nd x 1] vector of rotations (in degrees)
%       Svals : [nd x 1] vector of scaling factors
%       Kvals : [(nd*(nd-1)) x 1] vector of shear factors
% Output: 
%       A     : [nd+1 x nd+1] affine matrix

nd = 3; % only valid for 3D transformations
iK = [2,3,5,7,9,10];
if (length(Tvals)==nd) && (length(Rvals)==nd) && ...
                (length(Svals)==nd) && (length(Kvals)==length(iK)) && ...
                all(Svals>0)
    A = eye(4);
    % Shear
    A(iK) = Kvals;
    % Scale
    A(1:nd+2:nd*(nd+1)) = Svals;
    % Rotations:
    for i = 1:3
        if Rvals(i)~=0
            R = eye(nd);
            R(1:2,1:2) = [ cosd(Rvals(i)) , -sind(Rvals(i)) ; ...
                           sind(Rvals(i)) ,  cosd(Rvals(i))  ];
            if i>1
                R = circshift(R,[i,i]);
            end
            R(nd+1,nd+1) = 1;
            A = R*A;
        end
    end
    % Translation
    A(1:nd,end) = Tvals;
end
