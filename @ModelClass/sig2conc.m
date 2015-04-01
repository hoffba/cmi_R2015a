% ModelClass function
function sig2conc(self)
if ~self.gdchk && (length(self.ydata)>self.dceOpts.nbase)
    cf = cos(self.dceOpts.flipa);
    sf = sin(self.dceOpts.flipa);
    So = mean(self.ydata(1:self.dceOpts.nbase)) * ...
        (1 - exp(-self.dceOpts.TR * self.dceOpts.R1t) * cosd(self.dceOpts.flipa))...
        / ((1 - exp(-self.dceOpts.TR * self.dceOpts.R1t)) * sin(self.dceOpts.flipa));
    self.ydata = ((-1/self.dceOpts.TR) * log((1-(self.ydata./So)/sf) ...
                ./ (1-((self.ydata./So)/sf)*cf)) - self.dceOpts.R1t) / self.dceOpts.r1;
    self.gdchk = true;
end