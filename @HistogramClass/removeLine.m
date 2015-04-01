% HistogramClass function
% Remove histogram from plot
function removeLine(self,label)
% If label is not input, all histograms are cleared
if (nargin == 2)
if strcmp(label(end),'*') % Accepts wildcards at end of label
    label(end) = [];
    n = strncmp(label,self.labels,length(label));
else
    n = find(strcmp(label,self.labels),1);
end
else
n = 1:length(self.hLines);
end
delete(self.hLines(ishandle(self.hLines(n))));
self.hLines(n) = [];
self.binVals(:,n) = [];
self.labels(n) = [];
self.means(n) = [];
self.stdevs(n) = [];
if isempty(self.binVals)
self.setActive(false);
end