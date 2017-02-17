% CMIclass function
% Draw ROI on current slice
function drawROI(self,~,~)
% Draw ROI on current image
% roifun specifies the add/replace/remove function applied to the current mask
%       = 'Auto' -> determines whether to add/remove based on what's selected
%       = 'Replace'
%       = 'Add'
%       = 'Intersect'
%       = 'Remove'
if self.img.check && self.guicheck
    % Make sure there is an image showing
    if ~ishandle(self.haxes)
        self.dispUDview; % if not, rei-initialize the figure
    end
    % If overlay is showing, temporarily turn it off (can't have AlphaData)
    if (self.prmcheck || self.overcheck)
        overmod = true;
        tprmcheck = self.prmcheck;
        tovercheck = self.overcheck;
        self.prmcheck = false;
        self.overcheck = false;
        self.dispUDslice;
    else
        overmod = false;
    end
    % Draw ROI on slice
    if self.drawMode
        h = imfreehand(self.haxes);
        BW = h.createMask(self.hiover);
        h.delete;
    else
        axes(self.haxes);
        BW = roipoly;
    end
    % Return to overlay mode
    if overmod
        self.prmcheck = tprmcheck;
        self.overcheck = tovercheck;
    end
    if any(BW(:))
        self.img.mask.merge(self.roifun,BW,self.slc(self.orient),self.orient);
        self.dispUDmask;
        self.dispUDhist;
    end
end