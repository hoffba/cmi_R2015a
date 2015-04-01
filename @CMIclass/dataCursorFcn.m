% CMIclass function
% Data Cursor text function
function output_txt = dataCursorFcn(self,~,event_obj)
% function output_txt = myfunction(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
d = self.img.voxsz(1:3); d(self.orient) = [];
pos = fliplr(get(event_obj,'Position'));
matpos = ceil(pos ./ d);
pos(pos == 0) = 1;
dimg = self.img.getSlice(self.orient,self.vec,self.slc(self.orient));
val = dimg(matpos(1),matpos(2));
output_txt = {['X: ',num2str(pos(2),4)],...
              ['Y: ',num2str(pos(1),4)],...
              ['Z: ',num2str((self.slc(self.orient)-0.5)*self.img.voxsz(self.orient))],...
              ['matX: ',num2str(matpos(2),4)],...
              ['matY: ',num2str(matpos(1),4)],...
              ['Val: ',num2str(val)]};