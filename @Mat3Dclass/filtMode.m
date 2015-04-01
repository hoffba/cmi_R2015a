% Mat3Dclass function
function filtMode(self,wd)
% Filters the PRM non-zero values by a moving-mode function
%   to reduce random noise in the image
% Input: wd = 2-element vector (2D filtering only)

t = tic;
if nargin<2
    wd = [3,3];
end
str = 'Performing 2D Mode Filter on PRM ... ';
hw = waitbar(0,str);
for i = 1:self.dims(3)
    timg = self.mat(:,:,i);
    if any(timg(:))
        timg(timg==0) = nan;
        self.mat(:,:,i) = colfilt(timg,wd,'sliding',@mode);
    end
    waitbar(i/self.dims(3),hw,[str,num2str(i)]);
end
close(hw);
t = toc(t);
s = round(rem(t,60));
m = floor(s/60);
disp(['2D mode filter took ',num2str(m),':',num2str(s)]);