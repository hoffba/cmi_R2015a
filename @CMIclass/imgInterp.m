% CMIclass function
% Interpolate Image
function imgInterp(self,~,~)
if self.img.check
    d = self.img.dims(1:3);
    d0 = inputdlg('Interpolate to what dimensions?',...
                      'Interpolate',1,{num2str(d)});
    if ~isempty(d0)
        d0 = str2num(d0{1});
        if ~any(isnan(d0)) && (length(d0)==3)
            self.img.interp(d0(:));
            tdims = self.img.dims(1:3);
            ns = tdims(self.orient);
            self.slc = round(tdims/2);
            % Update GUI objects
            currslc = self.getSlc;
            set(self.h.slider_slice,'Value',currslc,...
                                    'SliderStep',[1,1]/(ns-1),...
                                    'Max',ns);
            set(self.h.edit_slc,'String',num2str(currslc));
            set(self.h.text_ns,'String',num2str(ns));
            self.dispUDslice;
        end
    end
end