function mask = cmi_kmask(r,type)
% Generate 3D kernel mask for filters
% Inputs:
%   r = window size radius (kernel = [2*r+1]
%   type = (string) type of window, e.g. 'circ'=3D-sphere

if isvector(r) && ischar(type) && ismember(type,{'circ','box'})
    N = zeros(1,3);
    N(1:length(r)) = r;
    switch type
        case 'box'
            mask = ones(2*N+1); % duh
        case 'circ'
            N = ceil((N-1)/2);
            [xx,yy,zz] = ndgrid(-N(1):N(1),-N(2):N(2),-N(3):N(3));
            N(N==0) = 1; % avoid divie by zero errors
            mask = sqrt(xx.^2/N(1)^2 + yy.^2/N(2)^2 + zz.^2/N(3)^2) <= 1;
    end
end

