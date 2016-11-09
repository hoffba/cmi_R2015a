% CMIclass function
% Create Montage of current image
function genMontage(self,x,~)
% Inputs:
%   
if self.img.check && self.dispcheck
    if isa(x,'matlab.ui.container.Menu')
    	F = self.grabFrames;
    elseif isstruct(x) && all(isfield(x,{'vdim','mdim','mind'}))
        F = self.grabFrames(x);
    else
        error('CMIclass/genMontage : Invalid inputs.')
    end
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