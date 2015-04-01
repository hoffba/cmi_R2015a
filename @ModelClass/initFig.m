% ModelClass function
function initFig(self)
if ~strcmp(self.fitOpts.Display,'off')
    % Define tags for plot objects:
    figstr = 'FitFig';
    lstr1 = 'FitLineData';
    lstr2 = 'FitLine';
    ftitle = 'FitTitle';
    newchk = isempty(self.h);
    % Create figure if needed
    if newchk || isempty(self.h.fig) || ...
            ~(ishandle(self.h.fig) && strcmp(get(self.h.fig,'Tag'),figstr))
        self.h.fig = figure('Tag',figstr);
        self.h.ax = axes('Parent',self.h.fig);
    end
    % Create data line if needed
    if newchk || isempty(self.h.dline) || ...
            ~(ishandle(self.h.dline) && strcmp(get(self.h.dline,'Tag'),lstr1))
        hold(self.h.ax,'on');
        self.h.dline = plot(self.h.ax,self.x,self.ydata,'.b','Tag',lstr1);
        hold(self.h.ax,'off');
    else
        set(self.h.dline,'YData',self.ydata);
    end
    % Create fit line if needed
    if newchk || isempty(self.h.fline) || ...
            ~(ishandle(self.h.fline) && strcmp(get(self.h.fline,'Tag'),lstr2))
        hold(self.h.ax,'on');
        self.h.fline = plot(self.h.ax,self.x,self.getDefYData,'-r','Tag',lstr2);
        hold(self.h.ax,'off');
    end
    % Create title with parameter values
    if newchk || isempty(self.h.title) || ...
            ~(ishandle(self.h.title) && strcmp(get(self.h.title,'Tag'),ftitle))
        self.h.title = title(self.h.ax,'');
        set(self.h.title,'Tag',ftitle);
    end
    set(self.h.title,'UserData',self.defs.(self.getModType)(self.cmodel(self.mtype)).labels')
end
