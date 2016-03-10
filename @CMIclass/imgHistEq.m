% CMIclass function
function imgHistEq(self,~,~)

stat = self.img.imgHistEq(self.vec);
if stat
    [imin,imax] = self.getColorMinMax;
    self.clim(self.vec,:) = [imin,imax];
    % update GUI objects:
    if self.guicheck
        [tmin,tmax] = self.getColorLimits;
        cpad = (imax - imin) / 2;
        set(self.h.slider_cmin,'Min',(imin-cpad),'Max',(imax+cpad),...
                                'Value',tmin,'Enable','On');
        set(self.h.slider_cmax,'Min',(imin-cpad),'Max',(imax+cpad),...
                                'Value',tmax,'Enable','On');
        set(self.h.edit_cmin,'String',num2str(tmin),'Enable','On');
        set(self.h.edit_cmax,'String',num2str(tmax),'Enable','On');
    end
    self.dispUDslice;
end