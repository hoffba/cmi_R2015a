% HistogramClass function
% Add vals to current histogram
function addVals(self,vals,label,varargin)
self.initPlotFig;
if ((nargin == 3) && ~isempty(self.binLocs)) % Make sure vals was input
    % Histogram the input values
    thistn = hist(vals(:),self.binLocs)';
    n = find(strcmp(label,self.labels),1);
    if isempty(n) % Add new histogram
        self.labels{end+1} = label;
        self.binVals(:,end+1) = thistn;
        self.hLines(end+1) = plot(self.binLocs,thistn/max(thistn),...
            'DisplayName',label,'Parent',self.haPlot);
        self.means(end+1) = mean(vals(:));
        self.stdevs(end+1) = std(vals(:));
        self.meds(end+1) = median(vals(:));
    else % Add vals to specified histogram
        self.binVals(:,n) = self.binVals(:,n) + thistn;
        set(self.hLines(n),'YData',self.binVals(:,n)/max(self.binVals(:,n)));
        oldn = sum(self.binVals(:,n));
        oldMean = self.means(n);
        newn = numel(vals);
        newMean = oldMean*oldn/(newn+oldn) +...
            mean(vals(:))*newn/(newn+oldn);
        % Update set statistics
        self.means(n) = newMean;
        self.stdevs(n) = sqrt((oldn*(self.stdevs(n)^2 + oldMean^2))...
                                / (oldn+newn) - newMean^2);
        tsum = 0;
        ind = 0;
        tot = sum(self.binVals(:,n));
        % ** Approximate median **
        while tsum<tot/2
            ind = ind+1;
            tsum = tsum+self.binVals(ind,n);
        end
        self.meds(n) = self.binLocs(ind);
    end
    self.dispUDstats;
    self.dispUDvisible;
end