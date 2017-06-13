% CMIclass function
% Clear mask
function clearMask(self,slc,sview)
% Inputs: slc - (not input) --> entire VOI
%             - 0 --> current slice in current view
%             - # in slice range
%         dim - slice dimension (if empty assumes current view)

if self.img.mask.check
    allchk = false;
    if (nargin==3) && ishandle(slc)
        tag = get(slc,'Tag');
        if strcmp(tag,'tools_voi_clearAll')
            allchk = true;
        elseif strcmp(tag,'tools_voi_clearSlice')
            sview = self.orient;
            slc = self.slc(sview);
        end
    elseif nargin==1
        allchk = true;
    end
    if (nargin<3) || isempty(sview)
        sview = self.orient;
    else
        sview = sview(1);
    end
    if (nargin<2) || isempty(slc)
        slc = self.slc(sview);
    else
        slc = slc(1);
    end
    if allchk
        self.img.mask.clear;
    else
        self.img.mask.clear(slc,sview);
    end
    self.dispUDmask;
end