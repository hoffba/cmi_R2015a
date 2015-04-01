function M = pts2T(phom,pref,Tflag,cor)
% Inputs:
%   X = homologous points [n x 3]
%   U = reference points [n x 3]
%   Tflag = type of transform
%           0: Translation
%           1: Rigid Body
%           2: Similarity
%           3: Affine

if (nargin<3)
    Tflag = 3;
end
if (nargin<4)
    cor = mean(phom,1);
end
M = eye(4);
npts = min(size(phom,1),size(pref,1));
if (nargin>1) && ismember(Tflag,0:3) && (size(phom,2)==3) && (size(pref,2)==3) ...
        && (size(phom,1)==size(pref,1))
    
    switch Tflag
        case 0
            nmin = 1;
        case 1
            nmin = 3;
        case 2
            nmin = 3;
        case 3
            nmin = 4;
    end
    if npts<nmin
        error('Not enough points for transform')
    end
    
    % Compute centroid-based coordinates
    phom = bsxfun(@minus,phom(1:npts,:),cor);
    pref = bsxfun(@minus,pref(1:npts,:),cor);
    
    if Tflag==3 % Full Affine
        
        % Linear least squares solution for full affine:
        pref = [pref,ones(npts,1)];
        phom = [phom,ones(npts,1)];
        M = [pref'*phom/(phom'*phom);0,0,0,1];

    else
        % Reference: 
        
        % Determine point-set centers
        chom = mean(phom,1);
        cref = mean(pref,1);
        
        if Tflag==0
            A = eye(3);
            T = (cref-chom)';
        else
            % Determine rotation matrix via SVD:
            X = bsxfun(@minus,pref,cref);
            U = bsxfun(@minus,phom,chom);
            [u,s,v] = svd(X' * U);
            Z = diag([1,1,det(u*v')]);
            R = v*Z*u';
            if Tflag==2
                % Determine scales
                S = trace(s*Z)/sum(sum(phom.^2))*(npts-1);
            else
                S = 1;
            end
            T = (cref' - (S*R)\chom');
            A = (S*eye(3)*R)';
        end
        M = [A,T;0,0,0,1];
    end
end