% CMIclass function
function scale2SNR(self,~,~)
% Scales image to SNR using current mask assumed to be noise

% Select which images in 4D to scale:
if self.img.dims(4)>1
    [ind,ok] = listdlg('ListString',self.img.labels);
    if ~ok
        return;
    end
else
    ind = 1;
end

N = mean(self.img.getMaskVals(ind))/1000;
self.img.imgScale(ind,1./N,zeros(1,length(N)),true);
self.clim(ind,:) = self.clim(ind,:)./repmat(N',1,2);
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