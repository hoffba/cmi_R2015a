% CMIclass function
function scale2HU(self,~)
% Scales image to HU using current mask
% Mask must contain:
%   2 separate regions over air and water
%   OR 1 region over water (assumes air is 0)

for i = 1:self.img.dims(4)
    om = self.img.scaleM(i);
    ob = self.img.scaleB(i);
    [m,b] = self.img.scale2HU(i);
    self.clim(i,:) = m/om * (self.clim(i,:)-ob) + b;
end
self.dispUDimg;

% update GUI objects:
if self.guicheck
    %[tmin,tmax] = self.getColorLimits;
    [imin,imax] = self.getColorMinMax;
    cpad = (imax - imin) / 2;
    set(self.h.slider_cmin,'Min',(imin-cpad),'Max',(imax+cpad),...
                            'Value',self.clim(self.vec,1),'Enable','On');
    set(self.h.slider_cmax,'Min',(imin-cpad),'Max',(imax+cpad),...
                            'Value',self.clim(self.vec,2),'Enable','On');
    set(self.h.edit_cmin,'String',num2str(self.clim(self.vec,1)),'Enable','On');
    set(self.h.edit_cmax,'String',num2str(self.clim(self.vec,2)),'Enable','On');
end