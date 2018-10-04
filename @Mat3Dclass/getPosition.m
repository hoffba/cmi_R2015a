% Mat3Dclass function
function p = getPosition(self,str,ind)

p=[];
if ~ischar(str)
    error('Invalid input.');
end

switch lower(str)
    case 'origin'
        p = self.slcpos;
    case 'direction'
        p = self.dircos;
    case 'coordinates'
        if nargin<2
            ind = [];
        elseif ~isnumeric(ind) || size(ind,2)~=3
            error('Invalid input.');
        end
        n = cross(self.dircos(1:3),self.dircos(4:6));
        M = [reshape(self.dircos,3,2),n',self.slcpos';zeros(1,3),1];
        d = self.dims(1:3);
        np = size(ind,1);
        if isempty(ind)
            np = prod(d);
            [X,Y,Z] = meshgrid((0:d(2)-1),(0:d(1)-1),(0:d(3)-1));
            ind = reshape(cat(4,X,Y,Z),np,3);
            clear X Y Z
        end
        p = M * [ (ind*diag(self.voxsz([2,1,3])))';ones(1,np)];
        p = p(1:3,:)';
end