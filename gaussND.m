function H = gaussND(s,r)
% Creates an N-D Gaussian filter
% Inputs: s = standard deviation (voxels, determines # dimensions)
%         r = (optional) radius of filter kernel

H = [];
if (nargin>0) && isnumeric(s)
    R = zeros(1,3);
    nd = length(s);
    if (nargin==1) || ~isnumeric(r)
        r = ceil((ceil(6*s-1)-1)/2);
    end
    if length(r)~=nd
        r = ceil(r(1))*ones(1,nd);
    end
    dim = cell(1,nd);
    for i = 1:3
        if i<=nd
            dim{i} = -r(i):r(i);
        else
            dim{i} = 0;
        end
    end
    [X{1},X{2},X{3}] = ndgrid(dim{:});
    R(1:nd) = r;
    H = ones(2*R+1);
    for i = 1:3
        if i<=nd
            H = H .* exp( -X{i}.^2 / (2*s(i)^2) ) / sqrt(2*pi*s(i));
        end
    end
    H = H/sum(H(:));
end