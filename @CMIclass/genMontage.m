% CMIclass function
% Create Montage of current image
function hf = genMontage(self,x,nrows)
% Inputs:
%   
if self.img.check && self.guicheck
    if isa(x,'matlab.ui.container.Menu')
    	F = self.grabFrames;
        nrows = [];
    elseif isstruct(x) && all(isfield(x,{'vdim','mdim','mind'}))
        F = self.grabFrames(x);
        if nargin<3
            nrows = [];
        end
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
        if isempty(nrows)
            nrows = round(str2double(inputdlg('Desired number of rows:','Rows',1,...
                {num2str(round(sqrt(size(timg,4))))})));
        end
        hf = figure;
        if nf>1
            montage(timg/255,'Size',[nrows,NaN]);
        else
            imshow(timg/255);
        end
    end
end