function X = filtGaussSep(X,str)
% Performs Gaussian filters over all dimensions of ND matrix X
% Input: X = matrix of values
%        str = smoothing strength
if (nargin==2) && isnumeric(X) && isnumeric(str)
    n = max(5,2*str);
    if n<3
        n = 3;
    end
    pn = round(n/2);
    X = padarray(X,pn*ones(1,3),'replicate');
    d = size(X);
    H = gaussdesign(1/str,2*n);
    for i = 1:length(d)
        if d(1)>50
            X = reshape(conv(X(:),H,'same'),d);
        end
        X = shiftdim(X,1);
        d = size(X);
    end
    X = X(pn+1:end-pn,pn+1:end-pn,pn+1:end-pn);
end
