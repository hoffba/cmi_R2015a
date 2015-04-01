% CMIclass function
% Load image
function stat = loadImg(self,x,fname)
% appcheck = self.appendcheck;

% Validate inputs:
if (nargin==3) && ishandle(x)
    switch get(x,'Tag')
        case 'file_openimg'
            appcheck = false;
        case 'file_append'
            appcheck = true;
    end
elseif (nargin>=2) && ~ishandle(x) && (isnumeric(x) || islogical(x))
    appcheck = logical(x); % input overrides property
else
    appcheck = self.appendcheck;
end
if (nargin<3) || ~(iscellstr(fname) || ischar(fname))
    fname = {};
end
tvec = self.img.dims(4); % need in case of append
stat = self.img.loadImg(appcheck,fname);
if stat % check that new image was actually loaded
    if appcheck
        d4 = self.img.dims(4);
        if self.applyallcheck
            vmin = ones(d4-tvec,1) * self.clim(self.vec,1);
            vmax = ones(d4-tvec,1) * self.clim(self.vec,2);
        else
            [vmin,vmax] = self.img.getColorMinMax(tvec+1:d4);
        end
        self.clim(tvec+1:d4,:) = [vmin vmax];
        self.vec = tvec + 1;
    else
        self.vec = 1;
        self.bgvec = 1;
        self.orient = 3;
        self.slc = round(self.img.dims(1:3)/2);
        [vmin,vmax] = self.img.getColorMinMax;
        self.clim = [vmin vmax];
        % Set GUI view buttons
        if self.guicheck
            self.h.button_slc.Value = 1;
            set([self.h.button_col,self.h.button_row],'Value',0);
        end
        self.activatePRM(false);
        notify(self,'LoadImgEvent');
    end
    % Update GUI objects:
    self.GUIupdate;
end