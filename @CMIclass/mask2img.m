% CMIclass function
function mask2img(self,str,~)

if self.img.check
    stat = self.img.mask2img(str);
    if stat
        self.clim(end+1,:) = [0,1];
        set(self.h.slider_4D,'Max',self.img.dims(4),...
                             'SliderStep',[1,1]/(self.img.dims(4)-1),...
                             'Enable','on');
        set(self.h.edit_vec,'Enable','on');
        set(self.h.text_nv,'String',num2str(self.img.dims(4)));
        self.dispUDslice;
    end
end