% CMIclass function
% Data Cursor text function
function output_txt = dataCursorFcn(self,~,event_obj)
% function output_txt = myfunction(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

% d = self.img.dims(1:3);
i = 1:3;
i(self.orient) = [];
pos = round(get(event_obj,'Position') + 0.5);
ind([i,self.orient]) = [ pos(2) , pos(1) , self.slc(self.orient) ];
vpos = self.img.getImageCoords(ind);
dimg = self.img.getSlice(self.orient,self.vec,self.slc(self.orient))';

val = dimg(pos(1),pos(2));
output_txt = {['X: ',num2str(vpos(1),4)],...
              ['Y: ',num2str(vpos(2),4)],...
              ['Z: ',num2str(vpos(3),4)],...
              ['matX: ',num2str(ind(1),4)],...
              ['matY: ',num2str(ind(2),4)],...
              ['matZ: ',num2str(ind(3),4)],...
              ['Val: ',num2str(val)]};