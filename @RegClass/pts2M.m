% RegClass function
function M = pts2M(self,varargin)
% M = pts2M;
% M = pts2M(refpts,hompts);
% M = pts2M(refpts,hompts,Ttype)

M = [];
if (nargin==1) || isa(varargin{1},'matlab.ui.control.UIControl')
    refpts = [self.points{1},ones(size(self.points{1},1),1)]*self.cmiObj(1).img.orient';
    hompts = [self.points{2},ones(size(self.points{2},1),1)]*self.cmiObj(2).img.orient';
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
% if (length(unique(refpts(:,3)))==1)
%     % 2D --> 3D for mapping photo to 3D image
%     M = pts2T(hompts(:,1:2),refpts(:,1:2),val);
%     A = [ M(1:2,1:2) , zeros(2,1) ; 0 , 0 , 1 ] ;
%     T = [ M(1:2,3) ; mean(hompts(:,3))-mean(refpts(:,3)) ];
%     M = [ A , T ; 0 0 0 1 ];
% else
    M = pts2T(hompts(:,[2,1,3]),refpts(:,[2,1,3]),val); % ref/hom backwards because Elastix works on inverse
    A = M(1:3,1:3);
    T = M(1:3,4);
% end

tpars = [reshape(A',1,[]),T'];
self.elxObj.setTx0(tpars,self.cmiObj(1).img.voxsz([2,1,3]),...
                             self.cmiObj(1).img.dims([2,1,3]),...
                             self.cmiObj(1).img.orient([2,1,3,4],[2,1,3,4]),...
                             'DefaultPixelValue',self.T0defVal);

self.showTx0;
self.setTchk(true);

