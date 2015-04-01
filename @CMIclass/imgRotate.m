% CMIclass function
% Rotate Image
function imgRotate(self,~,~)
if self.img.check
    answer = inputdlg('Rotate image CLOCKWISE by (degrees)');
    a = str2num(char(answer));
    if ~isempty(a) && isnumeric(a)
        a = rem(a,360);
        n = round(a/90);
        a = a - n*90;
        ncheck = mod(a,90);
        if (ncheck==0) % Quick rotation +/-90 or 180 degrees
            dorder = self.img.rotate90(self.orient,-n);
            self.slc = self.slc(dorder);
        else % non-90 rotation
            od = self.img.dims(1:3);
            self.img.rotate(self.orient,-a);% Default is CounterClockwise, so use negative value
            nd = self.img.dims(1:3);
            self.slc = self.slc(1:3) + round((nd-od)/2);
        end
        self.dispUDview;
    end
end