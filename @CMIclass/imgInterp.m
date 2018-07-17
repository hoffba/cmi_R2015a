% CMIclass function
% Interpolate Image
function imgInterp(self,d0,~)
if self.img.check
    if isa(d0,'matlab.ui.container.Menu')
        answer = inputdlg('Interpolate to what dimensions?',...
                          'Interpolate',1,{num2str(self.img.dims(1:3))});
        d0 = str2num(answer{1});
    end
    if isnumeric(d0) && (length(d0)==3)
        if ~any(isnan(d0)) && (length(d0)==3)
            self.img.interp(d0(:));
            tdims = self.img.dims(1:3);
            ns = tdims(self.orient);
            self.slc = round(tdims/2);
            % Update GUI objects
            currslc = self.getSlc;
            if ns>1
                set(self.h.slider_slice,'Value',currslc,...
                                        'SliderStep',[1,1]/(ns-1),...
                                        'Max',ns);
            else
                set(self.h.slider_slice,'Enable','off');
            end
            set(self.h.edit_slc,'String',num2str(currslc));
            set(self.h.text_ns,'String',num2str(ns));
            self.dispUDslice;
        end
    end
end