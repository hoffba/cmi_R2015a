% CMIclass function
% Calculate ADC map
function calcADC(self)
stat = self.img.adcCalc;
if stat
    self.clim(3,:) = [0 2500];
    self.vec = self.img.dims(4);
    self.dispUDslice;
end