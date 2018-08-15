% CMIclass function
% Rotate Image
function imgRotate(self,varargin)
if self.img.check
    if nargin==1 || isa(varargin{1},'matlab.ui.container.Menu')
        answer = inputdlg('Rotate image CLOCKWISE by (degrees)');
        a = str2double(char(answer));
        if isnan(a)
            error('Invalid input.')
        end
    elseif isscalar(varargin{1})
        a = varargin{1};
    else
        error('Invalid input.');
    end
    a = rem(a,360);
    n = round(a/90);
    a = a - n*90;
    if n~=0 % Quick rotation +/-90 or 180 degrees
        dorder = self.img.rotate90(self.orient,-n);
        self.slc = self.slc(dorder);
    end
    if a~=0 % non-90 rotation
        od = self.img.dims(1:3);
        self.img.rotate(self.orient,a);% Default is CounterClockwise, so use negative value
        nd = self.img.dims(1:3);
        self.slc = self.slc(1:3) + round((nd-od)/2);
    end
    self.dispUDview;
end