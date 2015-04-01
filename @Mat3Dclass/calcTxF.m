% ImageClass function
% Calculate transformation matrix
function T = calcTxF(self,U,V)
% U       = [n x 3] spatial coordinates of original space
% X       = [n x 3] spatial coordinates of transformed space
npts = length(U);
if self.check && (nargin==3) && (npts==length(V))
    % performs an affine transform on 3D image data
    U = [U ones(npts,1)]';
    V = [V ones(npts,1)]';
    if rank(U)>=4
        T = V*U'/(U*U');
    %     disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    %     disp(num2str(T))
    %     disp(['det: ' num2str(det(T))])
    %     disp(['norm: ' num2str(norm(T))])
    %     disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    else
        error('At least 4 non-collinear points needed to infer full affine transform.  Please select new points.');
    end
end