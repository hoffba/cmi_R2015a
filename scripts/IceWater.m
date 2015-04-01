% CMI script
function IceWater(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings

if cmiObj.img.dims(4) == 3
    bval = [120,1200];
    adcChk = true;
    ni = 1;
    scale = .001;
elseif cmiObj.img.dims(4) == 4
    bval = [3.8715115953,1108.5046458,1108.5046458,1195.33808874];
    adcChk = false;
    ni = 3;
    scale = 1000;
end
ns = cmiObj.img.dims(3);
out = zeros(ni*(ns+1),2);

for j = 1:ni
    ii = (j-1)*(ns+1);
    if adcChk
        adc = cmiObj.img.mat(:,:,:,3);
    else
        adc = log(cmiObj.img.mat(:,:,:,1)./cmiObj.img.mat(:,:,:,j+1))/(bval(j+1)-bval(1));
    end
    for i = 1:ns
        timg = adc(:,:,i);
        vals = timg(cmiObj.img.mask.mat(:,:,i));
        out(i+ii,:) = [mean(vals),iqr(vals)];
    end
    vals = adc(cmiObj.img.mask.mat);
    out(ii+ns+1,:) = [mean(vals),iqr(vals)];
end

assignin('base','out',out*scale)