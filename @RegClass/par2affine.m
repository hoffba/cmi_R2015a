% RegClass function
function M = par2affine(self,x)
% Convert sub-affine parameters to affine matrix

M = eye(4);
n = length(x);
switch n
    case 3 % Translation
        M(1:3,4) = x(:);
    case {6,7} % Euler
        if n==7 % Similarity
            s = x(7);
        else
            s = 1;
        end
        sf = sin(x(1:3));
        cf = cos(x(1:3));
        M(1:3,1:3) = [1,0,0;0,cf(1),-sf(1);0,sf(1),cf(1)] ...
                    *[cf(2),0,sf(2);0,1,0;-sf(2),0,cf(2)] ...
                    *[cf(3),-sf(3),0;sf(3),cf(3),0;0,0,1] ...
                    *(s*eye(3));
        M(1:3,4) = x(4:6);
    case 12 % Affine
        M = [reshape(x(1:9),3,3)',x(10:12)';0,0,0,1];
    otherwise
        error('RegClass:par2affine - Invalid number of parameters input.');
end
self.showTx0;