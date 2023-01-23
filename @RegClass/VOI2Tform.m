% RegClass function
function M = VOI2Tform(self,~,~)
% Use Ref/Hom VOIs to infer initial transform
%   * Scale and Translation only
M = [];
if self.cmiObj(1).img.mask.check && self.cmiObj(2).img.mask.check
        
    % Find spatial coordinates of mask voxels
    sz = zeros(2,3);
    cent = zeros(2,3);
    for i = 1:2
        d = self.cmiObj(i).img.dims(1:3);
        [vx,vy,vz] = ind2sub(d,find(self.cmiObj(i).img.mask.mat));
        n = length(vy);
        xyz = [vx,vy,vz,ones(n,1)] * self.cmiObj(i).img.orient';
        xyz(:,4) = [];
        sz(i,:) =  max(xyz) - min(xyz); % extent
        cent(i,:) = sum(xyz)/n; % centroid
    end
    S = sz(2,:)./sz(1,:);
    T = cent(2,:) - S.*cent(1,:);
    x = [ reshape(diag(S),1,[]) , T ];
%     x = [ reshape(diag(sz(2,:)./sz(1,:)),1,[]) , diff(cent) ];

    if strcmp(self.elxObj.outfmt,'.mhd')
        x = x([2,1,3,4],[2,1,3,4]);
    end
    
% OLD VERSION - from MHD use of image-centric coordinates assumed
%     % Scale / Translate based on VOI limits
%     voxsz1 = self.cmiObj(1).img.voxsz;
%     voxsz2 = self.cmiObj(2).img.voxsz;
%     d1 = self.cmiObj(1).img.dims(1:3);
%     d2 = self.cmiObj(2).img.dims(1:3);
%     [v0y,v0x,v0z] = ind2sub(d1,find(self.cmiObj(1).img.mask.mat));
%     [v1y,v1x,v1z] = ind2sub(d2,find(self.cmiObj(2).img.mask.mat));
%     v0x = ([min(v0x),max(v0x)] - d1(1)/2) * voxsz1(1);
%     v0y = ([min(v0y),max(v0y)] - d1(2)/2) * voxsz1(2);
%     v0z = ([min(v0z),max(v0z)] - d1(3)/2) * voxsz1(3);
%     v1x = ([min(v1x),max(v1x)] - d2(1)/2) * voxsz2(1);
%     v1y = ([min(v1y),max(v1y)] - d2(2)/2) * voxsz2(2);
%     v1z = ([min(v1z),max(v1z)] - d2(3)/2) * voxsz2(3);
%     x([1,5,9]) = [ diff(v1x)/diff(v0x) ,...
%                    diff(v1y)/diff(v0y) ,...
%                    diff(v1z)/diff(v0z) ];
%     x(isnan(x)) = 1; % correct for divide by zero
%     x(10:12) =   [ (sum(v1x)-sum(v0x)*x(1))/2 ,...
%                    (sum(v1y)-sum(v0y)*x(5))/2 ,...
%                    (sum(v1z)-sum(v0z)*x(9))/2 ];
               
    self.elxObj.setTx0(x,self.cmiObj(1).img.voxsz,self.cmiObj(1).img.dims(1:3),...
        self.cmiObj(1).img.orient,'DefaultPixelValue',self.T0defVal);
    self.setTchk(true);
    
    if self.guicheck
        set(self.h.popup_Transforms,'Value',5); % Set to Affine
        self.showTx0;
    end
    M = [ reshape(x,3,4) ; 0 0 0 1 ];
end
