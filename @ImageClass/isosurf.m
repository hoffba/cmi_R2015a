% ImageClass function
% Displays isosurface of either image data or mask
function isosurf(self,vec,val)
% vec = vector # of desired image (0 = mask)
if self.check && nargin==3 && isnumeric(vec) && isnumeric(val) && vec>=0
    vec = round(vec);
    ok = true;
    if vec
        fv = isosurface(self.mat(:,:,:,vec),val);
    elseif self.mask.check
        fv = isosurface(self.mask.mat,0.2);
    else % mask surface selected but no mask exists
        ok = false;
    end
    if ok
        figure;
        p = patch(fv,'EdgeColor','none','FaceColor',[228 209 192]/255);
        if vec
            isonormals(self.mat(:,:,:,vec),p);
        else
            isonormals(self.mask.mat,p);
        end
        daspect(1./self.voxsz);
        camlight
        view(122,16); axis tight off
    end
end