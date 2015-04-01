%ImageClass function
function plot4D(self)
if self.check && self.mask.check
    vals = self.getMaskVals;
    vals = mean(vals);
    figure
    plot(1:self.dims(4), vals)
end
