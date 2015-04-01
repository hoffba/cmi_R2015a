% Generate 2D histogram and Mutual Information
function [H,I] = hist2(x,y,xloc,yloc)
% Inputs:   x/y         = vectors of equal length
%           xloc/yloc   = bin locations for x and y
% Outputs:  H = 2D histogram bin values
%           I = Mutual Information metric

H = []; I = [];
if (nargin>1)
    nx = numel(x);
    if (numel(y)==nx)
        if (nargin<3)
            xnbin = 64;
            xloc = linspace(min(x),max(x),xnbin);
        else
            xnbin = length(xloc);
        end
        if (nargin<4)
            ynbin = 64;
            yloc = linspace(min(y),max(y),ynbin);
        else
            ynbin = length(yloc);
        end
        [xn, xbin] = histc(x,xloc);
        [yn, ybin] = histc(y,yloc);
        ind = ((xbin>0) & (ybin>0));
        Hbin = (xbin(ind)-1)*ynbin + ybin(ind);
        H = reshape(histc(Hbin,0.5:(xnbin*ynbin)),[xnbin,ynbin]);
        if nargout==2
            [xn,yn] = meshgrid(xn/nx,yn/nx);
            H = H/nx;
            I = H .* log(H ./ xn ./ yn);
            I = sum(I(~isnan(I)));
        end
    end
end


