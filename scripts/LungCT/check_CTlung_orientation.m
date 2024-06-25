function TF = check_CTlung_orientation(img,voxsz)
% Find orientation of shoulder bones to see if permute is needed

% BAH 2023-02-20                BW = max(img(i).mat(:,:,round(img(i).info.d(3)*4/5):end)>250,[],3);

% Grab top 1/4 slab to catch shoulders:
d = size(img);
BW = img(:,:,round(d(3)*3/4):end);

% Threshold to top 5% of values over 0HU
se = strel('sphere',3);
BW = max(imclose(BW > prctile(BW(BW>0),95),se),[],3);

% Find orientation
[y,x] = find(BW);
xy = [y*voxsz(2),x*voxsz(1)];
C = cov(xy);
[V,D] = eig(C);
D = diag(D);
[~,ind] = max(D);
TF = abs(V(ind,1)) > abs(V(ind,2));
