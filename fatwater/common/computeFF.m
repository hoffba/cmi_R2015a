function ff = computeFF( w,f,noiseBiasCorr )

if nargin<3
    noiseBiasCorr = 0;
end

denom = (abs(f) + abs(w));
denom2 = denom;
denom2(denom==0) = 1; % To avoid divide-by-zero issues
ff = 100*abs(f)./denom2;

if noiseBiasCorr>0
    fatregions = ff>50;
    watregions = ff<=50;
    denom2 = abs(f + w);
    denom2(denom==0) = 1; % To avoid divide-by-zero issues
    ff(watregions) = 100 - 100*abs(w(watregions))./denom2(watregions);
    ff(fatregions) = 100*abs(f(fatregions))./denom2(fatregions);
end