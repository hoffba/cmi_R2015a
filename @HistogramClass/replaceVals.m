% HistogramClass function
% Replace a histogram
function replaceVals(self,oldLabel,vals,newLabel)
if (nargin > 2)
    self.initPlotFig;
    if (nargin < 4)
        newLabel = oldLabel;
    end
    if strcmp(oldLabel(end),'*')
        oldLabel = oldLabel(1:end-1); % Remove the '*'
        n = find(strncmp(oldLabel,self.labels,(length(oldLabel)-1)));
    else
        n = find(strcmp(oldLabel,self.labels),1);
    end
    if isempty(n)
        self.addVals(vals,newLabel);
    else
        n = n(1);
        thistn = hist(vals(:),self.binLocs)';
        self.binVals(:,n) = thistn;
        set(self.hLines(n),'DisplayName',newLabel,...
            'YData',thistn(2:end-1)/max(thistn(2:end-1)),...
            'XData',self.binLocs(2:end-1));
        self.means(n) = mean(vals(:));
        self.stdevs(n) = std(vals(:));
        self.meds(n) = median(vals(:));
    end
    self.dispUDstats;
    self.dispUDvisible;
end