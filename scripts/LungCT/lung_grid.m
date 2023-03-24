function grid_label = lung_grid(mask,voxsz)

numpatch = 10000;
min_patch_sz = [3,3,3];

d = size(mask);
grid_label = uint16(zeros(d));

% Find mask limits:
[row,col,slc] = ind2sub(d,find(mask));
rmin = min(row);
rmax = max(row);
cmin = min(col);
cmax = max(col);
smin = min(slc);
smax = max(slc);
d_sub = [rmax-rmin,cmax-cmin,smax-smin]+1;
fov_sub = d_sub.*voxsz;

% Determine patch size
sz = (prod(fov_sub)/numpatch)^(1/3);
szd = ceil(sz*ones(1,3)./voxsz);
szd = max(szd,min_patch_sz);
N = ceil(d_sub./szd);

% Build the grid
C = num2cell(reshape(1:prod(N),N));
C = cellfun(@(x)x*ones(szd),C,'UniformOutput',false);
C = cell2mat(C);

% Insert grid into mask sized label matrix
d_C = N.*szd;
grid_label((1:d_C(1))+rmin-1,(1:d_C(2))+cmin-1,(1:d_C(3))+smin-1) = C;