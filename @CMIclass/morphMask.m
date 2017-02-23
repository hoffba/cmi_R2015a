% CMIclass function
% Perform morphological operation on mask (dilate/erode/open/close)
% morphMask(mtype,R)
% Inputs:
%   mtype = character describing what type of morphological operation:
%           d = dilate
%           e = erode
%           o = open (erode then dilate)
%           c = close (dilate then erode)
%   R     = radius of dilation [xy,z]
function morphMask(self,varargin)
if self.img.mask.check % First make sure there's something to morph
    % Check for manual/GUI inputs
    if (nargin==1) || isa(varargin{1},'matlab.ui.container.Menu')
        prompt = {'Operator ([d]ilate/[e]rode/[o]pen/[c]lose):',...
                  'Dilation Radius (pix):',...
                  'Dilation Height (0:2D, >0:3D):'};
        def = {'d','3','0'}; % 3D dilation, x-y radius = 3, z radius = 1
        answer = inputdlg(prompt,'Morph VOI:',1,def);
        if isempty(answer)
            return;
        else
            mtype = validatestring(answer{1},{'dilate','erode','open','close'});
            xyR = str2double(answer{2});
            zR = str2double(answer{3});
            if isnan(xyR) || isnan(zR)
                error('Invalid morph radius.')
            end
            R = [xyR,zR];
        end
    elseif nargin==3 && ischar(varargin{1}) && isvector(varargin{2})
        mtype = validatestring(varargin{1},{'dilate','erode','open','close'});
        R = varargin{2};
    else
        error('Invalid inputs.')
    end
    self.img.mask.morph(mtype,R);
    if self.guicheck
        self.dispUDmask;
    end
end