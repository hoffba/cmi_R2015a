% CMIclass function
% Data Cursor text function
function output_txt = dataCursorFcn(self,~,event_obj)
% function output_txt = myfunction(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

voxsz = self.img.voxsz(1:3); voxsz(self.orient) = [];

% Grab axes rectilinear XY coordinates
pos = get(event_obj,'Position');

% Convert to matrix indices
matpos = flip(ceil(pos ./ voxsz));
pos(pos == 0) = 1;

% Grab image value
dimg = self.img.getSlice(self.orient,self.vec,self.slc(self.orient));
val = dimg(matpos(1),matpos(2));

output_txt = {['X: ',num2str(pos(1),4)],...
              ['Y: ',num2str(pos(2),4)],...
              ['Z: ',num2str((self.slc(self.orient)-0.5)*self.img.voxsz(self.orient))],...
              ['matX: ',num2str(matpos(1),4)],...
              ['matY: ',num2str(matpos(2),4)],...
              ['Val: ',num2str(val)]};