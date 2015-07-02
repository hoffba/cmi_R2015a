% CMIcalss function
% Retrieve current slice
function varargout = getProp(self,varargin)
% Fetch parameter using current CMIclass settings:
% i.e. orient-related properties
% ** Number of return values should equal number of requested values.

varargout = cell(1,nargout);
% General settings:
a = self.orient;
if self.bgcheck
    v = self.bgvec;
else
    v = self.vec;
end
% Loop over properties to return:
for i = 1:nargout
    switch varargin{i}
        case 'fov'
            tval = self.img.dims(1:3).*self.img.voxsz;
        case 'slc'
            tval = self.slc(a);
        case 'nSlices'
            tval = self.img.dims(a);
        case 'nVec'
            tval = self.img.dims(4);
        case 'VOImajAxis'
            tval = self.img.getMaskMajorAxis(a);
        case 'ValLim'
            [tval(1),tval(2)] = self.img.getColorMinMax(v);
        case 'cLim'
            tval = self.clim(v,:);
        case 'cMap'
            if self.bgcheck
                tval = self.bgcmap;
            else
                tval = self.cmap;
            end
        case 'vStats'
            [tval(1),tval(2),tval(3),tval(4),tval(5)] = self.img.getStats(v);
        otherwise
            tval = [];
    end
    varargout{i} = tval;
end
if (nargin-1)~=nargout
    warning('CMIclass::getProp - Number of outputs not equal to requested values.')
end
