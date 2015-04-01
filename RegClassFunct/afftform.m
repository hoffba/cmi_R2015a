% Affine Transform Methods
function new_img = afftform(old_img,m,interpm,varargin)
% U = handles.ctrlpts{1}; % control points on reference
% X = handles.ctrlpts{2}; % control points on homologous
if size(varargin,2) == 1
    T = varargin{1};
elseif size(varargin,2) == 2
    U = varargin{1};
    X = varargin{2};
else
    errordlg('No inputs to transform')
end

% meshgrid method
if m==1
    if ~exist('T','var')
        npts = min([size(X,1) size(U,1)]);
        X = [X(1:npts,:) ones(npts,1)]';
        U = [U(1:npts,:) ones(npts,1)]';
        if rank(X)>=4
            disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            T = U*X'/(X*X');
            disp(num2str(T))
            disp(['det: ' num2str(det(T))])
            disp(['norm: ' num2str(norm(T))])
            disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        else
            error('At least 4 non-collinear points needed to infer full affine transform.  Please select new points.');
        end
    end
    tic
    new_img = im3affine(old_img,T,interpm);
    toc
    
    % affine_transform C-coded method
elseif m==2
    if ~exist('T','var')
        U = [U(:,2)-128 U(:,1)-128 U(:,3)-size(old_img,3)/2];
        X = [X(:,2)-128 X(:,1)-128 X(:,3)-size(old_img,3)/2];
        npts = min([size(X,1) size(U,1)]);
        X = [X(1:npts,:) ones(npts,1)]';
        U = [U(1:npts,:) ones(npts,1)]';
        if rank(X)>=4
            disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            T = inv(U*X'/(X*X'));
            disp(num2str(T))
            disp(['det: ' num2str(det(T))])
            disp(['norm: ' num2str(norm(T))])
            disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        else
            error('At least 4 non-collinear points needed to infer full affine transform.  Please select new points.');
        end
    end
    new_img = affine_transform(old_img,T);
    
    
    % affine w/ C-coded interpolation
elseif m==3
    if ~exist('T','var')
        U = [U(:,2)-128 U(:,1)-128 U(:,3)-size(old_img,3)/2];
        X = [X(:,2)-128 X(:,1)-128 X(:,3)-size(old_img,3)/2];
        npts = min([size(X,1) size(U,1)]);
        X = [X(1:npts,:) ones(npts,1)]';
        U = [U(1:npts,:) ones(npts,1)]';
        if rank(X)>=4
            disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            T = inv(U*X'/(X*X'));
            disp(num2str(T))
            disp(['det: ' num2str(det(T))])
            disp(['norm: ' num2str(norm(T))])
            disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        else
            error('At least 4 non-collinear points needed to infer full affine transform.  Please select new points.');
        end
    end
    new_img = affine(old_img,T);
else
    errordlg('not a transform method')
end