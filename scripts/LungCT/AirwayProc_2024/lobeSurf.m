function L_surf = lobeSurf(seg,voxsz)
% Generate surface triangulations for lung lobes

if ischar(seg)
    % Load segmentation from file
    if isfile(seg)
        [seg,~,fov,~,~,~] = cmi_load(1,[],seg);
        d = size(seg);
        voxsz = fov./d;
    else
        fprintf('File not found: %s\n',seg);
        return
    end
else
    d = size(seg);
    if nargin<2
        voxsz = ones(1,3);
    end
end

% Determine lobe codes
lobe = getLobeTags(seg);
lobe = lobe(cellfun(@(x)isscalar(x),{lobe.val})); % Use singular regions

nl = numel(lobe);
L_surf = cell(1,nl);
for i = 1:nl

    [xx,yy,zz] = ind2sub(d,find(seg==lobe(i).val));
    xyz = [xx,yy,zz].*voxsz; clear xx yy zz
    
    s = size(xyz);
    xyz = datasample(xyz,floor(s(1)/10)); % random sample of 10%
    xyz_b = boundary(xyz); % obtain indices for boundary set of voxels
    TR = triangulation(xyz_b,double(xyz));
    L_surf{i} = TR;
end
