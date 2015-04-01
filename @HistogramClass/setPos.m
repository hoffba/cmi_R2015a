% HistogramClass function
% Set histogram plot position
function setPos(self,pos)
if (nargin == 2)
    if ishandle(pos)
        self.hImgFig = pos;
        str = get(pos,'Units');
        if ~strcmp(str,'pixels')
            set(pos,'Units','pixels');
        end
        tpos = get(pos,'Position');
        set(pos,'Units',str);
        if (tpos(2) < 20) % Don't place plot under image if too low
            self.plotFigPos = get(0,'defaultfigureposition');
        else
            % Default histogram plot figure height = 200 (pixels)
            self.plotFigPos = [tpos(1) tpos(2)-280 tpos(3) 200];
        end
    elseif (length(pos)==4) && isnumeric(pos)
        self.plotFigPos = pos;
    end
end