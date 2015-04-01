% CMIclass function
function loc = getPoint(self,n,~)
% Retrieve location of user-selected point on image

loc = [];
if self.img.check && ishandle(self.hfig)
    if (nargin~=2) || isempty(n) || ~isnumeric(n)
        n = str2double(inputdlg('How many points?','# Pts',1,{'1'}));
    end
    if ~(isempty(n) || isnan(n))
        z = (self.slc(self.orient)-0.5)*self.img.voxsz(self.orient);
        figure(self.hfig);
        [x,y] = ginput(n(1));
        loc = [x,y,z*ones(n(1),1)];
    end
end