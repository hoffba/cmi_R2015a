function img = threshWT(img,wttype,lvls,thresh,opt)
% Performs 2D image filtering using thresholding in wavelet domain
% Syntax:
%   img = anidiff(img,);
% Inputs:
    
% Wavelete transform:
[C,S] = wavedec2(img,lvls,wttype);

% Thresholding of detail images:
i = S(1,1)*S(2,2)+1;
tC = C(i:end);
ii = tC>thresh;
tC(~ii) = 0;
if opt
    tC(ii) = sign(tC(ii)).*(abs(tC(ii))-thresh);
end
C(i:end) = tC;

% Wavelet reconstruction:
img = waverec2(C,S,wttype);
