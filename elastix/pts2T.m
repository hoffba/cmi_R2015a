function M = pts2T(U,X,Tflag)
% phom = M * pref
% Inputs:
%   X = homologous points [n x ndim]
%   U = reference points [n x ndim]
%   Tflag = type of transform
%           1: Translation
%           2: Rigid Body
%           3: Similarity
%           4: Affine

if (nargin<3)
    Tflag = 4;
end
ndim = size(U,2);
M = eye(ndim+1);
npts = min(size(U,1),size(X,1));

if (nargin>1) && ismember(Tflag,1:4) && (size(U,1)==size(X,1))
    
    % Make sure using same number of points each:
    U = U(1:npts,:);
    X = X(1:npts,:);
    
    if Tflag==4 % Full Affine
                
        % Linear least squares solution for full affine:
        X = [X,ones(npts,1)];
        M = [U'*X/(X'*X);zeros(1,ndim),1];
        
    else
        % Reference: 
        
        % Determine point-set centers
        ybar = mean(U,1);
        xbar = mean(X,1);
        
        if Tflag==1
            A = eye(ndim);
            T = (xbar-ybar)';
        else
            
            % Determine rotation matrix via SVD:
            X = bsxfun(@minus,X,xbar);
            U = bsxfun(@minus,U,ybar);
            [u,s,v] = svd(U' * X);
            Z = diag([ones(1,ndim-1),det(u*v')]);
            R = u*Z*v';
            
            if Tflag==3
                % Determine scales
                S = trace(s*Z)/sum(sum(X.^2)) * eye(ndim);
            else
                S = 1;
            end
            
            % Determine translation:
            A = S*R';
            T = xbar' - A*ybar';
        end
        M = [A,T;zeros(1,ndim),1];
    end
end