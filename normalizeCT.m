function Inorm = normalizeCT(I,varargin)
% Normalization function based on Gallardo-Estrella et al. 2016
% Inorm = normalizeCT(I)
% Inorm = normalizeCT(I,voxsz)
% Inorm = normalizeCT(I,BW)
% Inorm = normalizeCT(I,voxsz,BW)

BW = true;
n = 6;
if nargin>1
    if isscalar(varargin{1}) && isnumeric(varargin{1})
        n = varargin{1};
    elseif islogical(varargin{1}) && all(size(I)==size(varargin{1}))
        BW = varargin{1};
        if (nargin==3) && isscalar(varargin{2}) && isnumeric(varargin{2})
            n = varargin{2};
        end
    else
        error('Invalid input.');
    end
end

sig = [0,2.^(0:(n-2))];
d = size(I);
Inorm = I;
for i = 1:n % Loop over kernels
    
    Inorm = Inorm + lam*F;
end
