% ImageClass function
function [m,b] = scale2HU(self,vec)
% Scales image to HU using current mask
% Mask must contain either:
%       - 2 separate regions over air and water
%       - or only a region over water

m = 1; b = 0;
tissue = [-1000,10]; % Water=0HU, soft tissue=10-40HU, air=-1000HU, fat=-120
if tissue(2) == 10
    'Using muscle (10HU)'
end
if self.check && self.mask.check
    if nargin==1
        vec = 1;
    end
    cc = bwconncomp(self.mask.mat);
    if cc.NumObjects < 3
        ntot = prod(self.dims(1:3));
        vals = [0,0];
        for i = 1:cc.NumObjects
            vals(i) = mean(self.mat(cc.PixelIdxList{1} + ntot*(vec-1)));
        end
        vals = sort(vals);
        om = self.scaleM(vec);
        ob = self.scaleB(vec);
        m = diff(tissue)/diff(vals)*om;
        b = m * (ob - vals(2));
        self.imgScale(vec,m,b);
    else
        error('Input VOI is not valid. Need two regions, air and water.');
    end
end