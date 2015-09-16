function c = calcConfluence(I)
% c = calcConfluence(I)
% Calculates global confluence of binary image I

hw = waitbar(0,'Blurring image');

% Blur the binary image to get density map:
I = filtGaussSep(double(I),15);

% Calculate Euler at multiple density thresholds
dth = 0.05;
th = 0:dth:1;
nth = length(th);
X = zeros(1,nth);
% waitbar(0,hw,'Calculating Euler over density thresholds ...');
for i = 1:nth
    X(i) = imEuler3dEstimate(I>th(i));
%     waitbar(i/nth,hw);
end
% delete(hw);
% figure,plot(th,X)

% Calculate confluence (clumpiness)
X = (X-1).^2;
c = sqrt( dth/2 * (X(1) + X(end) + 2*sum(X(2:end-1)) ) );
