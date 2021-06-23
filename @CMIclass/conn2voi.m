% CMIclass function
function conn2voi(self,h,~)
% Finds connected region in thresholded image to save to mask
%   ** uses the VOI Mode to replace/add/intersect/remove

tmask = [];
if (nargin==1) || (ishandle(h) && strcmp(get(h,'Tag'),'tools_connThresh')) ...
        || (ischar(h) && strcmpi(h,'thresh'))
    % Use image threshold as mask
    if self.chk2D
        tmask = self.img.getSlice(self.orient,self.vec,self.slc(self.orient));
    else
        tmask = self.img.mat(:,:,:,self.vec);
    end
    tmask = (tmask >= self.img.thresh(self.vec,1)) ...
          & (tmask <= self.img.thresh(self.vec,2));
elseif self.img.mask.check
    % Use current VOI as mask
    if self.chk2D
        tmask = self.img.mask.getSlice(self.orient,self.vec,self.slc(self.orient));
    else
        tmask = self.img.mask.mat;
    end
end

if ~isempty(tmask)
    % User input for point selection:
    set(self.hfig,'WindowStyle','modal');
    title(self.haxes,'Select Region:','Visible','on');
    % axes(self.haxes);
    [x,y] = getpts(self.haxes);
    title(self.haxes,'','Visible','off');
    set(self.hfig,'WindowStyle','normal');

    t = tic;

    % Convert 2D point to matrix index:
    order = 1:3; order(self.orient) = [];
    voxsz = self.img.voxsz; voxsz(self.orient) = [];
    mx = round(x/voxsz(2)+0.5);
    my = round(y/voxsz(1)+0.5);
    if self.chk2D
        d = size(tmask);
        ii = sub2ind(d,my(:),mx(:));
    else
        subs(:,[order,self.orient]) = [my(:),mx(:),ones(length(mx),1)*self.slc(self.orient)];
        ii = sub2ind(self.img.dims,subs(:,1),subs(:,2),subs(:,3));
    end
    ii(~tmask(ii)) = [];
    if isempty(ii)
        error('Selected points must be within threshold limits.')
    end

    % Find connected regions:
    hw = waitbar(0,'Finding connected regions ...');
    cc = bwconncomp(tmask);
    waitbar(0.5,hw,'Finding selected region ...');
    ct = 1;
    tmask(:) = false;
    while ~isempty(ii)
        ind = ismember(ii,cc.PixelIdxList{ct});
        if any(ind)
            tmask(cc.PixelIdxList{ct}) = true;
            ii(ind) = [];
        end
        ct = ct+1;
    end
    delete(hw);

    % Set VOI:
    if strcmp(self.roifun,'Auto')
        if self.chk2D
            omask = self.img.mask.getSlice(self.orient,1,self.slc(self.orient));
        else
            omask = self.img.mask.mat;
        end
        if self.img.mask.check && any(omask(:) & tmask(:))
            str = 'remove';
        else
            str = 'add';
        end
    else
        str = lower(self.roifun);
    end
    if self.chk2D
        self.img.mask.merge(str,tmask,self.slc(self.orient),self.orient);
    else
        self.img.mask.merge(str,tmask);
    end
    self.dispUDroi;

    toc(t)
end
