% RegClass function
function M = pts2M(self,varargin)
% M = pts2M;
% M = pts2M(refpts,hompts);
% M = pts2M(refpts,hompts,Ttype)

M = [];
if (nargin==1) || isa(varargin{1},'matlab.ui.control.UIControl')
    refpts = self.points{1};
    hompts = self.points{2};
    Ttype = get(self.h.popup_Transforms,'String');
    Ttype = Ttype{get(self.h.popup_Transforms,'Value')};
elseif (nargin<3)
    error('Invalid inputs. Need both reference and homologous points.');
elseif (size(x,2)~=3) || (size(y,2)~=3)
    error('Input point matrices need to be [N x 3].');
elseif (nargin<4)
    Ttype = 'Affine';
end

% Determine MIN number of points:
switch Ttype
    case 'Translation'
        val = 0;
        n = 1;
    case 'Euler'
        val = 1;
        n = 3;
    case 'Similarity'
        val = 2;
        n = 3;
    case 'Affine'
        val = 3;
        n = 4;
    case 'NONE'
        return;
    otherwise
        error(['Invalid Transform Type: ',Ttype]);
end

% Number of points must be the same
np = min(size(refpts,1),size(hompts,1));
if np<n
    error(['Number of points (',num2str(np),') less than MIN (',num2str(n),').']);
end
% Points must be centered for calculation (RegClass default to center of image)
off1 = self.cmiObj(1).img.voxsz .* self.cmiObj(1).img.dims(1:3) /2;
% off2 = self.cmiObj(2).img.voxsz .* self.cmiObj(2).img.dims(1:3) /2;
% refpts = refpts(1:np,:) - repmat(off1,np,1);
% hompts = hompts(1:np,:) - repmat(off2,np,1);
M = pts2T(refpts,hompts,val,off1); % ref/hom backwards because Elastix works on inverse
% M(1:3,4) = M(1:3,4) + off2';
A = M(1:3,1:3)';
T = M(1:3,4)';

self.elxObj.setTx0([A(:)',T],self.cmiObj(1).img.voxsz,...
                       self.cmiObj(1).img.dims(1:3),...
                       'DefaultPixelValue',self.T0defVal);

self.showTx0(M(1:3,:));
self.setTchk(true);

