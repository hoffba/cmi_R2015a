function stat = initMain(self)
%Function for initializing all the callbacks for the RegClass

stat = false;
% Check if figure needs to be loaded
if isempty(self.h) && self.guicheck
    % All GUI objects and their callbacks are set here:
    self.regFig;
end
% Update GUI callbacks
if ~isempty(self.h)
    
    % ~~~~~~~ FIGURE POSITIONING ~~~~~~~~~
    % (Linux) figure size adjustments (pixels):
        brdsz = 2; % Border size
        headsz = 46; % Includes header, filemenubar, top border
    % Retrieve screen size for positioning:
        set(0,'units','pixels');
        bpos = get(0,'screensize');
    % Calculate CMI GUI positions:
        h = self.cmiObj(1).h.mainFig;
        set(h,'units','pixels');
        CMIpos = get(h,'Position');
        p1 = [brdsz , bpos(4)-CMIpos(4)-headsz , CMIpos(3:4)];
        p2 = [brdsz , bpos(4)/2-CMIpos(4)-headsz , CMIpos(3:4)];
    % Calculate imgReg GUI position:
        set(self.h.regFig,'units','pixels');
        pReg = get(self.h.regFig,'Position');
        pReg(1:2) = [bpos(3)-pReg(3)-brdsz , bpos(4)-pReg(4)-headsz];
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % ~~~~~~~~~~~~~~~ Set figure properties:
    set(self.h.regFig,   'CloseRequestFcn', @self.GUIcloseFcn,...
                         'Position',pReg);
    set(self.cmiObj(1).h.mainFig,  'Color', [0.5 0 0],...
                                   'Name', 'Fixed Image',...
                                   'Units','pixels',...
                                   'Position',p1);
    set(self.cmiObj(2).h.mainFig,  'Color', [0 0 0.5],...
                                   'Name', 'Moving Image',...
                                   'Units','pixels',...
                                   'Position',p2);
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            
    % Set listeners for plotting points on CMIclass figures
    addlistener(self.cmiObj,'slc','PostSet',    @self.plotPts);
    addlistener(self.cmiObj,'orient','PostSet', @self.plotPts);
    addlistener(self.cmiObj,'hiover','PostSet', @self.setBDFcn);
    addlistener(self,       'points','PostSet', @self.listenPoints);
    addlistener(self.elxObj,'Tx0','PostSet',@self.listenT0);
    
    % Set listener for initializing after loading new data sets
    addlistener(self.cmiObj,'LoadImgEvent',@self.initImgData);
    
    stat = true;
end

