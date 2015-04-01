% CMIclass function
% Create Montage of current image
function genMontage(self,~,~)
if self.img.check && self.dispcheck
    F = self.grabFrames;
    if ~isempty(F)
        nf = length(F);
        timg = zeros(size(F(1).cdata,1),size(F(1).cdata,2),3,nf);
        for i = 1:nf
            [im,map] = frame2im(F(i));
            if ~isempty(map)
                im = ind2rgb(im,map);
            end
            timg(:,:,:,i) = double(im);
        end
        warning('off','Images:initSize:adjustingMag');
        dim = sqrt(size(timg,4));
        figure;montage(timg/255,'Size',[round(dim),ceil(dim)]);
    end
end