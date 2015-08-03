function M = pts2T(Y,X,Tflag)
% phom = M * pref
% Inputs:
%   X = homologous points [n x ndim]
%   U = reference points [n x ndim]
%   Tflag = type of transform
%           0: Translation
%           1: Rigid Body
%           2: Similarity
%           3: Affine

if (nargin<3)
    Tflag = 3;
end
ndim = size(Y,2);
M = eye(ndim+1);
npts = min(size(Y,1),size(X,1));

if (nargin>1) && ismember(Tflag,0:3) && (size(Y,1)==size(X,1))
    
%     switch Tflag
%         case 0
%             nmin = 1;
%         case 1
%             nmin = 3;
%         case 2
%             nmin = 3;
%         case 3
%             nmin = 4;
%     end
%     if npts<nmin
%         error('Not enough points for transform')
%     end
    
    % Make sure using same number of points each:
    Y = Y(1:npts,:);
    X = X(1:npts,:);
    
    if Tflag==3 % Full Affine
        
        % Linear least squares solution for full affine:
        X = [X,ones(npts,1)];
        Y = [Y,ones(npts,1)];
        M = [Y(:,1:ndim)'*X/(X'*X);zeros(1,ndim),1];
        M(1:ndim,end) = -M(1:ndim,end);

    else
        % Reference: 
        
        % Determine point-set centers
        ybar = mean(Y,1);
        xbar = mean(X,1);
        
        if Tflag==0
            A = eye(ndim);
            T = (xbar-ybar)';
        else
            
            % Determine rotation matrix via SVD:
            X = bsxfun(@minus,X,xbar);
            Y = bsxfun(@minus,Y,ybar);
            [u,s,v] = svd(Y' * X);
            Z = diag([ones(1,ndim-1),det(u*v')]);
            R = u*Z*v';
            
            if Tflag==2
                % Determine scales
                S = trace(s*Z)/sum(sum(X.^2)) * eye(ndim);
            else
                S = 1;
            end
            
            % Determine translation:
            T = (S*R)*xbar' - ybar';
            A = S*R;
        end
        M = [A,T;zeros(1,ndim),1];
    end
end