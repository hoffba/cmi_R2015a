% This function finds the optimal Rigid/Euclidean transform in 3D space
% It expects as input a Nx3 matrix of 3D points.
% It returns R, t
% Written by Nghia Ho

% expects row data
function [R,t] = rotateTransFunct(A, B)
    if nargin ~= 2
	    error('Missing parameters');
    end    
    assert(size(A, 1) == size(B,1) && size(A,2) == size(B,2) && size(A, 3) == size(B,3));

    centroid_A = mean(A);
    centroid_B = mean(B);

    N = size(A,1);

    H = (A - repmat(centroid_A, N, 1))' * (B - repmat(centroid_B, N, 1));

    [U,~,V] = svd(H);

    R = V*U';

    if det(R) < 0
        disp('Reflection detected');
        V(:,3) = V(:,3)*-1;
        R = V*U';
    end

    t = -R*centroid_A' + centroid_B';
    
    %personal added code so output is 4x4 transformation matrix
    f = eye(4);
    f(1:3, 1:3) = R;
    R = f;
    n = eye(4);
    n(1:3, 4) = t;
    t = n;
    
end

