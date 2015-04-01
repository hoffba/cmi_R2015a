% CMIclass function
function autoScale(self,~,~)
% Automatically scales CT images based on histogram
[tmin,tmax] = self.getColorLimits;
m0 = self.img.scaleM(self.vec);
b0 = self.img.scaleB(self.vec);
stat = self.img.autoScale(self.vec);
% update GUI objects:
if stat && self.guicheck
    m = self.img.scaleM(self.vec);
    b = self.img.scaleB(self.vec);
    tmin = (tmin - b0)*m/m0 + b;
    tmax = (tmax - b0)*m/m0 + b;
    self.clim(self.vec,:) = [tmin,tmax];
    [imin,imax] = self.getColorMinMax;
    cpad = (imax - imin) / 2;
    set(self.h.slider_cmin,'Min',(imin-cpad),'Max',(imax+cpad),...
                            'Value',tmin,'Enable','On');
    set(self.h.slider_cmax,'Min',(imin-cpad),'Max',(imax+cpad),...
                            'Value',tmax,'Enable','On');
    set(self.h.edit_cmin,'String',num2str(tmin),'Enable','On');
    set(self.h.edit_cmax,'String',num2str(tmax),'Enable','On');
end