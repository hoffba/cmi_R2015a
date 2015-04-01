% ImageClass function
function D = getMaskMajorAxis(self,vdim)
% Returns length of major axis of mask for specified slice
if nargin==1
    vdim = 3;
end
D = zeros(self.dims(3),1);
for i = 1:self.dims(3)
    tmask = self.mask.getSlice(vdim,1,i);
    if ~isempty(tmask) && (nnz(tmask)>0)
        tvoxsz = self.voxsz;
        tvoxsz(vdim) = [];
        p = regionprops(tmask,'MajorAxisLength','Orientation','Area');
        ind = find(p.Area==max(p.Area),1);
        d = p.MajorAxisLength(ind);
        a = p.Orientation(ind);
        D(i) = d * sqrt((sind(a)*tvoxsz(1))^2 ...
                      + (cosd(a)*tvoxsz(2))^2);
    end
end
D = max(D);