function chi = fitfunct(self,refimg,homimg,mask,sx,sy,sz,thx,thy,thz,trx,tray,trz,shx,shy,shz)
%OPTIMIZER Summary of this function goes here
%   Detailed explanation goes here

timg = self.cmiObj(2).img.mat; 
mi = mutualinformation(timg, self.cmiObj(1).img.mat);
chi = -mi;

end

