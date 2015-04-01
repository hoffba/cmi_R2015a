% CMI script
function traceADC(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings

b = 585.95;

d = cmiObj.img.dims(1:3);
adc = zeros(d);
for i = 1:3
    adc = adc + log(cmiObj.img.mat(:,:,:,1)./cmiObj.img.mat(:,:,:,1+i))/b;
end
adc(cmiObj.img.mat(:,:,:,1)<1000) = 0;
tmat = cmiObj.img.mat;
tmat(:,:,:,end) = adc*10^6;
cmiObj.img.setMat(tmat);
%assignin('base','adc',adc/3);