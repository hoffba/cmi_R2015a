% FitClass function
% Fit model to data
function pout = fitData(self,idefs)
pout = [];
if self.OptimChk
    if nargin<2
        idefs = self.defs(self.mind);
    end
    if self.mind && ~isempty(self.ydata) && ...
            (length(self.xdata)==length(self.ydata)) && ...
            all(isfield(idefs,{'func','par0','bounds'}))
        if strcmp(self.fitOpts.Display,'iter')
            self.initFig;
        end
        pout = lsqnonlin(idefs.func,...
                         idefs.par0,...
                         idefs.bounds(1,:),... % lower
                         idefs.bounds(2,:),... % upper
                         self.fitOpts);
        if strcmp(self.fitOpts.Display,'final')
            self.initFig;
            self.updatePlot('data');
            self.updatePlot('fit',feval(idefs.func,pout));
            pause(0.01);
        end
    end
end