% ImageClass function
% Calculate ADC map from 2 b-values
function stat = adcCalc(self)
stat = false;
if self.dims(4) == 2
    answer = inputdlg('Input b-values:','Calculate ADC',1,{'120 1200'});
    bvals = str2num(answer{1});
    adc = 10^6 * log(((self.mat(:,:,:,1)-self.scaleB(1))/self.scaleM(1)) ./ ...
        ((self.mat(:,:,:,2)-self.scaleB(2))/self.scaleM(2))) / (bvals(2) - bvals(1));
    adc(adc<realmin) = realmin;
    self.mat = cat(4,self.mat,adc);
    self.dims(4) = 3;
    self.labels{3} = 'ADC [e-6 mm^2/2]';
    self.scaleM(3) = 1;
    self.scaleB(3) = 0;
    self.thresh(3,:) = 10^6*[-1 1];
    stat = true;
end