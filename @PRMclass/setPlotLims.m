% PRMclass function
function setPlotLims(self,xlim,ylim)

if ishandle(self.hascatter) && (nargin==3) ...
        && (numel(xlim)==2) && (numel(ylim)==2)
    
    set(self.hascatter,'XLim',xlim,'YLim',ylim)
    
    % Re-set threshold lines
    h = allchild(self.hascatter);
    h = h(strncmp(get(h,'Tag'),'Line',4));
    xdata = cell2mat(get(h,'XData'));
    ydata = cell2mat(get(h,'YData'));
    delete(h);
    m = diff(ydata,1,2)./diff(xdata,1,2);
    b = ydata(:,1) - m.*xdata(:,1);
    infi = isinf(m);
    b(infi) = xdata(infi,1);
    plotlines(self.hascatter,m,b);
end