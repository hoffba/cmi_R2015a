% FitClass function
function initFig(self)
if ~strcmp(self.fitOpts.Display,'off')
    % Define tags for plot objects:
    figstr = 'FitFig';
    lstr1 = 'FitLineData';
    lstr2 = 'FitLine';
    ftitle = 'FitTitle';
    newchk = isempty(self.h) || ~all(ishandle(self.h));
    
    % Create figure if needed
    if newchk
        self.h.fig = figure('Tag',figstr);
        self.h.ax = axes('Parent',self.h.fig);
        
        hold(self.h.ax,'on');
        self.h.dline = plot(self.h.ax,self.xdata,self.ydata,'.b','Tag',lstr1);
        hold(self.h.ax,'off');
        
        hold(self.h.ax,'on');
        self.h.fline = plot(self.h.ax,self.xdata,self.yfit,'-r','Tag',lstr2);
        hold(self.h.ax,'off');
        
        self.h.title = title(self.h.ax,'');
        set(self.h.title,'Tag',ftitle);
        set(self.h.title,'UserData',self.defs(self.mind).labels')
    else
        set(self.h.dline,'YData',self.ydata,'XData',self.xdata);
        set(self.h.fline,'YData',self.yfit,'XData',self.xdata);
    end
    
end
