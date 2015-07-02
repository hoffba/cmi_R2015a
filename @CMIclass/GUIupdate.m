% CMIclass function
function GUIupdate(self)
% Updates CMI GUI objects with current properties/limits

if self.guicheck
    buttonvals = zeros(1,3);
    buttonvals(self.orient) = 1;
    set(self.h.button_row,'Value',buttonvals(1));
    set(self.h.button_col,'Value',buttonvals(2));
    set(self.h.button_slc,'Value',buttonvals(3));
    cslc = self.slc(self.orient);
    td = self.img.dims([self.orient,4]);
    enslc = 'Off';
    if td(1)>1
        enslc = 'On';
    end
    envec = 'Off';
    if td(2)>1
        envec = 'On';
    end
    set(self.h.edit_slc,'String',num2str(cslc),...
                        'Enable',enslc);
    set(self.h.slider_slice,'Min',1,...
                            'Max',max(td(1),2),...
                            'Value',cslc,'Enable',enslc,...
                            'SliderStep',[1 1]/max(td(1)-1,1));
    set(self.h.text_ns,'String',num2str(td(1)));
    set(self.h.edit_vec,'String',num2str(self.vec),...
                        'Enable',envec);
    set(self.h.slider_4D,'Min',1,...
                         'Max',max(td(2),2),...
                         'Value',self.vec,'Enable',envec,...
                         'SliderStep',[1 1]/max(td(2)-1,1));
    set(self.h.text_nv,'String',num2str(td(2)));
    set(self.h.edit_veclabel,'String',self.img.labels{self.vec},...
                             'Enable','On');
    [vext,clim] = self.getProp('ValLim','cLim');
    cpad = diff(vext) / 2;
    set(self.h.slider_cmin,'Min',(vext(1)-cpad),'Max',(vext(2)+cpad),...
                            'Value',clim(1),'Enable','On');
    set(self.h.slider_cmax,'Min',(vext(1)-cpad),'Max',(vext(2)+cpad),...
                            'Value',clim(2),'Enable','On');
    set(self.h.edit_cmin,'String',num2str(clim(1)),'Enable','On');
    set(self.h.edit_cmax,'String',num2str(clim(2)),'Enable','On');
    set([self.h.slider_thMin,self.h.slider_thMax],...
        'Min',vext(1),'Max',vext(2));
    self.dispUDview;
end