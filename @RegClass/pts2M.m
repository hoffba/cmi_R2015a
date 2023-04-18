% RegClass function
function M = pts2M(self,varargin)
% M = pts2M;                        - uses points from RegClass
% M = pts2M(refpts,hompts);         - assumes rectilinear coordinates and affine T
% M = pts2M(refpts,hompts,Name/Value pairs)

M = [];

Tforms = {'Translation', 1; ...
          'Euler',       3; ...
          'Similarity',  3; ...
          'Affine',      4 };

if (nargin==1) || isa(varargin{1},'matlab.ui.control.UIControl')
    coord_str = 'rect';
    refpts = self.points{1};
    hompts = self.points{2};
    Ttype = get(self.h.popup_Transforms,'String');
    Ttype = Ttype{get(self.h.popup_Transforms,'Value')};
else
    % Parse inputs
    p = inputParser;
    addRequired(p,'refpts',@(x)isnumeric(x)&&size(x,2)==3);
    addRequired(p,'hompts',@(x)isnumeric(x)&&size(x,2)==3);
    addParameter(p,'coord_str','rect',@(x)ismember(x,{'ind','rect','ref'}));
    addParameter(p,'Ttype','Affine',@(x)ismember(x,Tforms(:,1)));
% elseif (nargin<3)
%     error('Invalid inputs. Need both reference and homologous points.');
% elseif (size(x,2)~=3) || (size(y,2)~=3)
%     error('Input point matrices need to be [N x 3].');
% elseif (nargin<4)
%     Ttype = 'Affine';
end

% Determine MIN number of points:
val = find(strcmp(Tforms(:,1),Ttype),1);
n = Tforms{val,2};

% Make sure number of points is not less than required for transform type
np = min(size(refpts,1),size(hompts,1));
if np<n
    error(['Number of points (',num2str(np),') less than MIN (',num2str(n),').']);
end

% Convert points from rectilinear to reference coordinates

% Reference:
np = size(refpts,1); 
fvoxsz = self.cmiObj(1).img.voxsz;
fdims = self.cmiObj(1).img.dims(1:3);
ffov = fdims.*fvoxsz;
fT = self.cmiObj(1).img.orient;
refpts = refpts / diag(fvoxsz([2,1,3])); % Convert to ijk indices
% refpts(:,2) = fdims(2) - refpts(:,2) + 1;        % Flip y-dimension
refpts = [refpts,ones(np,1)] * fT';      % Transform to referenced space

% Homologous
np = size(hompts,1); 
mvoxsz = self.cmiObj(2).img.voxsz;
mdims = self.cmiObj(2).img.dims(1:3);
mfov = mdims.*mvoxsz;
mT = self.cmiObj(2).img.orient;
hompts = hompts / diag(mvoxsz([2,1,3]));  % Convert to ijk indices
% hompts(:,2) = mdims(2) - hompts(:,2) + 1;      % Flip y-dimension
hompts = [hompts,ones(np,1)] * mT';       % Transform to referenced space

% if (length(unique(refpts(:,3)))==1)
%     % 2D --> 3D for mapping photo to 3D image
%     M = pts2T(hompts(:,1:2),refpts(:,1:2),val);
%     A = [ M(1:2,1:2) , zeros(2,1) ; 0 , 0 , 1 ] ;
%     T = [ M(1:2,3) ; mean(hompts(:,3))-mean(refpts(:,3)) ];
%     M = [ A , T ; 0 0 0 1 ];
% else
    M = pts2T(refpts(1:np,[2,1,3]),hompts(1:np,[2,1,3]),val);
    % M = pts2T(hompts(1:np,1:3),refpts(1:np,1:3),val);
    A = M(1:3,1:3);
    T = M(1:3,4);
% end

tpars = [reshape(A',1,[]),T'];
self.elxObj.setTx0(tpars,fvoxsz([2,1,3]),fdims([2,1,3]),fT,'DefaultPixelValue',self.T0defVal);

self.showTx0;
self.setTchk(true);

