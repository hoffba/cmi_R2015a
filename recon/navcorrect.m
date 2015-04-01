function dwi = navcorrect(imdat,navdat)
% performs a correction with navigator echo and FFT
% returns a 256 x 256 image matrix

np2 = size(imdat,2);
[nv,nnav] = size(navdat);
nvout = max(np2,nv);
zf1 = nvout/np2; % zero fill factor for the np dimension
zf2 = nvout/nv; % zero fill factor for the nv dimension

% zf1=256/np2;
% zf2=256/nv;

dc = mean(mean(imdat(:,(end-4):end)));

% Equate nav and img dims
navdata = zeros(zf1*np2,zf2*nv);
imgdata = zeros(zf1*np2,zf2*nv);
n0 = round((zf1*np2 - nnav)/2)+1;
n00 = round(((zf1-1)*np2)/2)+1;
n000 = round((zf2-1)*nv/2)+1;

% Roughly center nav and img echos
navdata(n0:n0+nnav-1,n000:n000+nv-1) = (navdat-dc)'; % Size = zf*(nimgdat by nv)
imgdata(n00:n00+np2-1,n000:n000+nv-1) = (imdat-dc)'; % Size = zf*((nimgdat by nv)
z1 = navdata(:,1:zf2*nv);
y = mean(z1,2);
x = ones(1,zf2*nv);
yy = y * x;
ref = fftshift(fft(yy));     % Note, this already has zero-filled dims
fz1 = fftshift(fft(navdata)); % Note, this already has zero-filled dims
imgdat1 = fftshift(fft(imgdata)); % FFT wrt readout only, this is zero-filled

% Try mean nav phase error instead of 1D phase error on readout axis
% Let ref amplitude be the weighting function of mean phase error
dc = mean((ref .* conj(fz1) ./ (abs(fz1) + 0.000001)),1);% Norm by abs(fz1) so dc amp is constant
dc = dc ./ (abs(dc)+0.000001);
x = ones(zf1*np2,1);
dcfix = x * dc;
imgdatc = imgdat1 .* dcfix;
imgdat1 = imgdatc';
imgdat = fft(imgdat1); % FFT wrt phase enc only, this is zero-filled
imgdat1 = imgdat';     % Back to original orientation
for ii=1:(zf1*np2)
    xx = fftshift(imgdat1(ii,:));
    imgdat1(ii,:) = xx;
end
imgdat = fliplr(imgdat1);
imgdat=flipdim(imgdat,1);
dwi(:,:) = abs(imgdat);


