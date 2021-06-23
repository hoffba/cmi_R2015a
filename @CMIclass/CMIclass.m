classdef CMIclass < handle
    properties (SetObservable, SetAccess=private, GetAccess=public)
        
        % Constants associated with ImageClass
        img                     % ImageClass object (created in contructor)
        orient = 3;             % Dimension orientation (3 = slice)
        
        clim                    % Main image color limits
        slc                     % slice numbers to show [row col slc]
        vec = 1;                % 4D index to show
        bgvec = 1;              % 4D index for background image
        roifun = 'Auto';        % Function for new ROIs
        drawMode = 1;           % VOI drawing function (0=roipoly, 1=imfreehand)
        
        histObj                 % HistogramClass object (created in constructor)
        
        % Display colors
        cmap = 'gray';          % Colormap for image
        bgcmap = 'gray';        % Colormap for background
        ncolors = 128;          % number of colors for display
        voispec = 'om';         % Define VOI edge type/color
        voimarksz = 2;          % Define VOI marker size
        thspec = 'oc';            % Define Threshold edge type/color
        thmarksz = 2;           % Define Threshold edge marker size
        dispPos                 % Position of display figure
        dalpha = 1;             % Transparency of overlay image (0:1)
        checkerd = 20;          % Checkerboard size
        checkerM = [];          % Checkerboard logical matrix
        
        % Display handles
        h                       % Handles to main figure (including children)
        hfig                    % Handle to display figure
        haxes                   % Handle of image axes
        hibg                    % Handle of background image
        hiover                  % Handle of overlay image
        hvoi                    % Handle of VOI plot
        hthresh                 % Handle of threshold plot
        hdcm                    % Handle of DataCursorMode object
        haplot                  % Handle of extra axes for plotting various data
        
        % Display parameters
        guicheck = false;       % Check whether GUI is displayed
        applyallcheck = false;  % Check for adjust all image colors
        histcheck = false;      % Determines whether to display histogram plots
        overcheck = false;      % Check for overlay display
        bgcheck = false;        % Check for adjusting background image
        prmcheck = false;       % Check to display PRM overlay
%         appendcheck = false;    % Check for whether to append new images
        
        % Other options
        chk2D = false; % Connectivity check (2D or 3D)
    end
    events
        LoadImgEvent
        ChangeView
    end
    methods % Constructor/Destructor
        % Constructor with handle inputs for Umich3 listeners
        function self = CMIclass(dchk)
            self.img = ImageClass(self);
            if nargin
                if isempty(dchk)
                    self.guicheck = false;
                elseif isa(dchk,'CMIclass')
                    % link to existing GUI
                    self.h = dchk.h;
                    self.guicheck = dchk.guicheck;
                else
                    self.guicheck = logical(dchk);
                end
            else
                self.guicheck = true;
            end
            self.initMain;
            if self.guicheck
                self.histObj = HistogramClass(self.h.mainFig);
            end
        end
        % Destructor method
        function delete(self)
            % Clean-up for all related figures, GUI, and member classes
            % Close figures
            if self.guicheck
                hfs = [self.hfig,self.h.mainFig];
                delete(hfs(ishghandle(hfs)));
            end
            % Delete member objects
            if isa(self.histObj,'HistogramClass') && isvalid(self.histObj)
                delete(self.histObj)
            end
            if isa(self.img,'ImageClass') && isvalid(self.img)
                delete(self.img)
            end
        end
    end
    methods (Access = private)
        % re-initialize the object
        function initialize(self)
            self.img.initialize;
            %self.histObj.setActive(false);
            self.orient = 3;
            self.vec = 1;
            self.bgvec = 1;
            self.prmcheck = false;
            self.overcheck = false;
            self.applyallcheck = false;
            self.dispUDcmap;
            self.dispUDhist;
        end
        
        % Listener for display figure position changes
        function changeDispPos(self,~,~)
            if ishandle(self.hfig)
                self.dispPos = get(self.hfig,'Position');
            end
        end
        
        % Initialize display figure
        function dispFigs(self)
            if ~isempty(self.hfig) && ishghandle(self.hfig)
                figure(self.hfig) % Figure exists, so just set it to current figure
            else
                if isempty(self.dispPos) % Only happens the first time the image is loaded
                    % Default figure 'Position' = [xleft ybottom 560width 420height]
                    % Toolbar height = 20 pixels
                    if ishandle(self.h.mainFig)
                        tunits = get(self.h.mainFig,'Units');
                        if ~strcmp(tunits,'pixels') % Must measure in pixels for this operation
                            set(self.h.mainFig,'Units','pixels');
                        end
                        tpos = get(self.h.mainFig,'Position');
                        set(self.h.mainFig,'Units',tunits); % Re-set to MainFig's default units
                        self.dispPos = [tpos(1)+tpos(3)+4 tpos(2)+tpos(4)-450 560 420];
                    else % Not created in "cmi" GUI, so set to arbitrary position
                        tpos = get(0,'ScreenSize');
                        self.dispPos = [(tpos(3)-560)/2 (tpos(4)-420)/2 560 420];
                    end
                end
                % Create MATLAB graphics objects for dispaying images
                self.hfig = figure('Position',self.dispPos,'Color',get(self.h.mainFig,'Color'),...
                                   'CloseRequestFcn','');
                self.haxes = axes('Visible','Off','Parent',self.hfig,'YDir','reverse',...
                                  'Position',[0,0,1,1]);
                self.hibg = image('Parent',self.haxes,'CData',[]);
                self.hiover = image('Parent',self.haxes,'CData',[]);
                % Update data cursor text function
                self.hdcm = datacursormode(self.hfig);
                set(self.hdcm,'UpdateFcn',@self.dataCursorFcn,'Enable','Off');
                % Set up the figure's colormap
                set(self.hfig,'ColorMap',[eval([self.bgcmap '(' num2str(self.ncolors) ')']);...
                                          eval([self.cmap   '(' num2str(self.ncolors) ')'])]);
                % Set up listener for figure position
                addlistener(self.hfig,'LocationChanged',@self.changeDispPos);
                % Setup toolbar options
                strs = {'New Figure','Open File','Edit Plot','Rotate 3D',...
                    'Brush/Select Data','Link Plot','Insert Legend','Insert Colorbar',...
                    'Hide Plot Tools','Show Plot Tools'};
                for i = 1:length(strs)
                    th = findall(self.hfig,'ToolTipString',strs{i});
                    set(th,'Visible','off');
                end
                % Initialize histogram object if necessary
                self.histObj.setPos(self.hfig);
            end
        end
        % Update background image
        function dispUDbg(self)
            if self.img.check && self.guicheck
                if ~ishandle(self.hfig)
                    self.dispUDview;
                % only update if background overlay is activated
                elseif (self.prmcheck || self.overcheck || ~isempty(self.checkerM))
                    % values are arbitrarily displayed, scaled to colormap
                    timg = (self.getImgSlice('img',self.bgvec) - self.clim(self.bgvec,1)) / ...
                        diff(self.clim(self.bgvec,:)) * (self.ncolors - 1);
                    timg(timg > (self.ncolors - 1)) = self.ncolors - 1;
                    set(self.hibg,'CData',timg);
                end
            end
        end
        % Update foreground image
        function dispUDimg(self)
            if self.img.check && self.guicheck
                if ~ishandle(self.hfig)
                    self.dispUDview;
                else
                    % Determine mask for transparency
                    if (self.prmcheck && self.img.prm.check)
                        tprm = self.getImgSlice('prm');
                        adata = ~isnan(tprm);
                        tprm(~adata) = 0;
                        timg = self.ncolors + tprm;
                        adata = double(adata & (tprm>0)) * self.dalpha;
                    elseif ~isempty(self.checkerM)
                        timg = self.getImgSlice('img');
                        adata = self.checkerM;
                        if self.img.mask.check
                            adata = adata & self.getImgSlice('mask');
                        end
                        adata = adata * self.dalpha;
                    elseif (self.overcheck && self.img.mask.check)
                        timg = self.getImgSlice('img');
                        adata = self.getImgSlice('mask');
                        adata = double( adata ...
                            & (timg>=self.img.thresh(self.vec,1)) ...
                            & (timg<=self.img.thresh(self.vec,2)) ) * self.dalpha;
                    else
                        timg = self.getImgSlice('img');
                        adata = 1; % No transparency
                    end
                    % Determine image to show
                    if ~self.prmcheck
                        tclim = self.clim(self.vec,:);
                        timg = (timg - tclim(1)) ./ diff(tclim) .* (self.ncolors - 1);
                        timg(timg < 1) = 1;
                        timg = timg + self.ncolors;
                    end
                    % Update front image
                    set(self.hiover,'CData',timg,'AlphaData',adata);
                    % ~~~~~~~~~~~~~~~~~~~~~~~
                    % Bug fix - sometimes overlays don't show (looks gray), this resets the image somehow
                    %axis on
                    %axis off
                    % ~~~~~~~~~~~~~~~~~~~~~~~ end fix
                end
            end
        end
        % Update ROI
        function dispUDroi(self)
            if self.img.check && self.guicheck
                if ~ishandle(self.hfig)
                    self.dispUDview;
                else
                    % Find ROI outline
                    if ~(self.overcheck || self.prmcheck)
                        [voir,voic] = find(edge(self.getImgSlice('mask'),'Canny'));
                        tsz = self.img.voxsz; tsz(self.orient) = [];
                        voir = tsz(1) * voir - tsz(1)/2;
                        voic = tsz(2) * voic - tsz(2)/2;
                    else
                        voir = [];
                        voic = [];
                    end
                    % Check to see if the current handle points to valid object
                    handlecheck = ~(isempty(self.hvoi) || ~ishandle(self.hvoi));
                    % Determine how to update the ROI plot
                    if (isempty(voir) && handlecheck)
                        delete(self.hvoi);
                    elseif (~handlecheck && ~isempty(voir))
                        hold(self.haxes,'on');
                        self.hvoi = plot(self.haxes,voic,voir,self.voispec,'MarkerSize',self.voimarksz);
                        hold(self.haxes,'off');
                    elseif ~isempty(voir)
                        set(self.hvoi,'XData',voic,'YData',voir,...
                            'MarkerSize',self.voimarksz,...
                            'Marker',self.voispec(1),...
                            'MarkerEdgeColor',self.voispec(2),...
                            'MarkerFaceColor',self.voispec(2));
                    end
                end
            end
        end
        % Update threshold plot
        function dispUDthreshplot(self)
            if self.img.check && self.guicheck
                if ~ishandle(self.hfig)
                    self.dispUDview;
                else
                    if ~(self.overcheck || self.prmcheck)
                        tmask = self.getImgSlice('img');
                        tmask = (tmask>=self.img.thresh(self.vec,1)) ...
                              & (tmask<=self.img.thresh(self.vec,2));
                        if nnz(tmask)==numel(tmask)
                            voir = [];
                            voic = [];
                        else
                            [voir,voic] = find(edge(tmask,'Canny'));
                            tsz = self.img.voxsz; tsz(self.orient) = [];
                            voir = tsz(1) * voir - tsz(1)/2;
                            voic = tsz(2) * voic - tsz(2)/2;
                        end
                    else
                        voir = [];
                        voic = [];
                    end
                    % Check to see if the current handle points to valid object
                    handlecheck = ~(isempty(self.hthresh) || ~ishandle(self.hthresh));
                    % Determine how to update the ROI plot
                    if (isempty(voir) && handlecheck)
                        delete(self.hthresh);
                    elseif (~handlecheck && ~isempty(voir))
                        hold(self.haxes,'on');
                        self.hthresh = plot(self.haxes,voic,voir,self.thspec,'MarkerSize',self.thmarksz);
                        hold(self.haxes,'off');
                    elseif ~isempty(voir)
                        set(self.hthresh,'XData',voic,'YData',voir,...
                            'MarkerSize',self.thmarksz,...
                            'Marker',self.thspec(1),...
                            'MarkerEdgeColor',self.thspec(2),...
                            'MarkerFaceColor',self.thspec(2));
                    end
                end
            end
        end
        % Update display image slice
        function dispUDslice(self)
            self.dispUDimg;
            th = findall(self.hfig,'ToolTipString','Data Cursor');
            if (self.overcheck || self.prmcheck || ~isempty(self.checkerM))
                self.dispUDbg;
                set(th,'Visible','off');
            else
                set(th,'Visible','on');
            end
            self.dispUDroi;
            self.dispUDthreshplot;
        end
        % Update display view
        function dispUDview(self)
            if isa(self.img,'ImageClass') && self.img.check && self.guicheck
                self.dispFigs;
                set(self.hfig,'Name',self.img.name)
                self.dispUDslice;
                self.setChecker;
                dpix = self.img.voxsz(1:3); dpix(self.orient) = [];
                tdims = self.img.dims(1:3); tdims(self.orient) = [];
                fov = dpix .* tdims;
                set(self.haxes,'XLim',[0 fov(2)-dpix(2)/8],'YLim',[0 fov(1)-dpix(1)/8],...
                    'DataAspectRatio',[1 1 1]);
                set([self.hiover self.hibg],'XData',[dpix(2) fov(2)]-dpix(2)/2,...
                                            'YData',[dpix(1) fov(1)]-dpix(1)/2);
            end
        end
        % Update display colormap
        function dispUDcmap(self)
            if ishandle(self.hfig)
                tcmap = eval([self.bgcmap '(' num2str(self.ncolors) ')']);
                if self.prmcheck
                    tcmap = [tcmap; self.img.prm.cmap];
                else
                    tcmap = [tcmap; eval([self.cmap '(' num2str(self.ncolors) ')'])];
                end
                set(self.hfig,'ColorMap',tcmap);
            end
        end
        % Update displays associated with the mask
        function dispUDmask(self)
            if self.img.check
                if self.prmcheck % need to re-calc PRM map
                    self.img.calcPRM(self.vec);
                end
                self.dispUDslice;
                self.dispUDhist;
            end
        end
        
        % Update histogram plots
        function dispUDhist(self,opt)
            % opt = 0: Volume
            %       1: Slice
            %       []: Both
            if self.guicheck
                tcheck = self.histcheck && self.img.mask.check;
                if (tcheck ~= self.histObj.active)
                    self.histObj.setActive(tcheck);
                end
                if (self.histObj.active && self.img.mask.check)
                    tvec = self.vec;
                    tmin = self.clim(tvec,1);
                    tmax = self.clim(tvec,2);
                    self.histObj.initialize(tmin,tmax,prod(self.img.voxsz));
                    if (nargin == 1)
                        opt = [0 1]; % No input means update both
                    end
                    if any(opt == 0) % Volume
                        timg = self.img.mat(:,:,:,tvec);
                        timg = timg(self.img.getThreshMask(tvec));
                        self.histObj.replaceVals('Volume',timg);
                    end
                    if any(opt == 1) % Slice
                        tslc = self.slc(self.orient);
                        timg = self.img.getSlice(self.orient,self.vec,tslc);
                        timg = timg(self.img.getThreshMask(tvec,self.orient,tslc));
                        self.histObj.replaceVals('Slice',timg);
                    end
                end
                figure(self.h.mainFig);
            end
        end

    end
end