function sseg = yacta2sext(seg)

useg = unique(seg); useg(1)=[];
n = numel(useg);
if n==2
    Rval = 10;
    Lval = 20;
elseif n==6
    Rval = [10,20,30];
    Lval = [40,50,60];
else
    sseg = [];
    return
end

d = size(seg);
sseg = zeros(size(seg));
Z = repmat(permute(1:d(3),[3,1,2]),d(1),d(2),1);

% Find RL sextants
L = ismember(seg,Rval);
ind = find(any(L,[1,2]));
ind = linspace(ind(1)-.5,ind(end)+.5,4);
sseg(L & (Z>=ind(1)) & (Z<ind(2))) = 63;
sseg(L & (Z>=ind(2)) & (Z<ind(3))) = 62;
sseg(L & (Z>=ind(3)) & (Z<ind(4))) = 61;

% Find RL sextants
L = ismember(seg,Lval);
ind = find(any(L,[1,2]));
ind = linspace(ind(1)-.5,ind(end)+.5,4);
sseg(L & (Z>=ind(1)) & (Z<ind(2))) = 66;
sseg(L & (Z>=ind(2)) & (Z<ind(3))) = 65;
sseg(L & (Z>=ind(3)) & (Z<ind(4))) = 64;