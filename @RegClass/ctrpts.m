% RegClass function
function p = ctrpts(self,ind)
% Converts points vector to image-centered points used in Elastix

p= [];
if (nargin==1)
    ind = 1:2; % Both
elseif ~isnumeric(ind) || any(~ismember(ind,1:2));
    error('RegClass:ctrpts - input must be 1 (Ref) or 2 (Hom), or both');
end
if ~isempty(ind)
    np = min(cellfun(@(x)size(x,1),self.points(ind)));
    p = zeros(np,3,length(ind));
    for i = 1:length(ind)
        ii = ind(i);
        toff = self.cmiObj(ii).img.voxsz.*self.cmiObj(ii).img.dims(1:3)/2;
        p(:,:,i) = self.points{ii}(1:np,:) - repmat(toff,np,1);
    end
end
