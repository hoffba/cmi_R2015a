function [L_surfs,lobe_ids] = getLobeSurfs(seg,voxsz)
% Finds distinct lung lobes and returns surface renderings for visualization

dim = size(seg);
lobe_ids = unique(seg(seg>0));
isoval = movmean([0;lobe_ids],2); isoval(1) = [];

nl = numel(lobe_ids);
L_surfs = cell(1,nl);
[X,Y,Z] = meshgrid((1:dim(1))*voxsz(1),(1:dim(2))*voxsz(2),(1:dim(3))*voxsz(3));
for i = 1:nl
    fv = isosurface(X,Y,Z,permute(seg,[2,1,3]),isoval(i));
    fv = reducepatch(fv,0.1);
    L_surfs{i} = triangulation(fv.faces,fv.vertices);
end












return



warning('off','MATLAB:triangulation:PtsNotInTriWarnId')

dim = size(seg);
lobe_ids = unique(seg);

nl = numel(lobe_ids);
L_surfs = cell(1,nl);
for i = 1:nl
    [xx,yy,zz] = ind2sub(dim,find(seg==lobe_ids(i)));
    tmp = [xx,yy,zz].*voxsz; clear xx yy zz
    
    s = size(tmp);
    tmp2 = datasample(tmp,floor(s(1)/10));
    tmp3 = boundary(tmp2);
    TR = triangulation(tmp3,double(tmp2));
    L_surfs{i} = TR;
end
