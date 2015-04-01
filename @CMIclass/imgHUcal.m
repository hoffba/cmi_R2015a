% CMIclass function
function imgHUcal(self,~,~)
% Calibrates and image to HU based on water and air in image

if self.img.check
    self.setView(1);
    
    tm = self.img.scaleM(self.vec);
    tb = self.img.scaleB(self.vec);
    
    % First select water region
    title(self.haxes,'Select WATER region:')
    axis on; axis off;
    axes(self.haxes)
    BW = roipoly;
    timg = self.img.getSlice(self.orient,self.vec,self.slc(self.orient));
    watADU = (mean(timg(BW)) - tb) / tm;
    
    % Next select air
    title(self.haxes,'Select AIR region:')
    axes(self.haxes)
    BW = roipoly;
    timg = self.img.getSlice(self.orient,self.vec,self.slc(self.orient));
    airADU = (mean(timg(BW)) - tb) / tm;
    title(self.haxes,'')
    
    % Scale image using these mean values
    % HU scaled up by 1000 (for fld)
    disp(['Air = ',num2str(airADU),' , Water = ',num2str(watADU)])
    m = 1000/(watADU - airADU);
    b = -airADU * m;
    self.img.imgScale(self.vec,m,b);
end