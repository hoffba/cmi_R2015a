% Generate affine matrix (2D/3D) from rotation, translation, scale, and shear factors
% M = makeAffTform(r,t,s,k)
% Inputs:   r = Rotation - [1xnD]
%           t = Translation - [1xnD]
%           s = Scale - [1xnD] / scalar for isotropic / [] for none
%           k = Shear - [1x(nD^2-nD)] / [] for none
% Output:   M = Affine matrix - [(nD+1) x (nD+1)]
function M = makeAffTform(r,t,s,k)

if (nargin>1) && (length(r)==length(t)) && any(length(r)==[2,3])
    d = length(t);
    chk3D = (d==3);
    if ~exist('s','var') || isempty(s)
        s = [];
    end
    if ~exist('k','var')
        k = [];
    end
    
    % Rotation
    cr = cosd(r);
    sr = sind(r);
    if chk3D
        M = [1 0 0 0 ; 0 cr(1) -sr(1) 0 ; 0 sr(1) cr(1) 0 ; 0 0 0 1] * ...
            [cr(2) 0 sr(2) 0 ; 0 1 0 0 ; -sr(2) 0 cr(2) 0 ; 0 0 0 1] * ...
            [cr(3) -sr(3) 0 0 ; sr(3) cr(3) 0 0 ; 0 0 1 0 ; 0 0 0 1];
    else
        M = [cr sr 0 ; -sr cr 0 ; 0 0 1];
    end

    % Scale
    if ~isempty(s)
        if length(s)==1
            s = s*ones(1,d);
        end
        M = M * diag([s,1]);
    end

    % Shear
    if length(k)==(d^2-d)
        if chk3D
            M = M * [1 k(1:2) 0 ; k(3) 1 k(4) 0 ; k(5:6) 1 0 ; 0 0 0 1];
        else
            M = M * [1 k(1) 0 ; k(2) 1 0 ; 0 0 1];
        end
    end

    % Translate
    M(1:d,end) = t;
end