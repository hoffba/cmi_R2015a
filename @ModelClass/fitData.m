% ModelClass function
% Fit model to data
function pout = fitData(self)
pout = [];
if self.OptimChk
    ok = true;
    if strcmp(self.getModType,'perf')
        self.calcSo; % Perfusion fits are done on S values
        if strcmp(self.getModName,'refreg') && (length(self.RRdata)~=length(self.ydata))
            ok = false; % RR ydata is not valid, cannot continue with fit
            disp('Need to set reference region values first!')
        end
    end
    if ok && ~isempty(self.ydata) && (length(self.x)==length(self.ydata))
        if strcmp(self.fitOpts.Display,'iter')
            self.initFig;
        end
        tdefs = self.getDefs;
        pout = lsqnonlin(@(p)self.fitFunc(p),...
                            tdefs.par0,...
                            tdefs.bounds(1,:),... % lower
                            tdefs.bounds(2,:),... % upper
                            self.fitOpts);
        if strcmp(tdefs.name,'BiExp') && (pout(3)>pout(2))
            % diffusion value order matters - higher first
            %   and f needs to be the slow fraction
            pout = pout([1,3,2,4]);
            pout(4) = 1 - pout(4);
        end
        if strcmp(self.fitOpts.Display,'final')
            self.initFig;
            self.updatePlot('data');
            self.updatePlot('fit',feval(@(x)self.(self.getDefs.func)(x),pout));
            pause(0.01);
        end
    end
end