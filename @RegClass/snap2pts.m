% RegClass function
function stat = oneshot(self,x,y,Ttype)
% stat = snap2pts;
%           Uses points already input into RegClass.points
% stat = snap2pts(refpts,hompts);
%           Uses refpts [n x 3] and hompts [n x 3]
% Use user-input points to infer Affine-type transform 
%       and apply to moving image

if (nargin==1) || ishandle(x)
    varargin = {};
elseif (nargin<3)
    error('Invalid inputs. Need both reference and homologous points.');
elseif (size(x,2)~=3) || (size(y,2)~=3)
    error('Input point matrices need to be [N x 3].');
else
    varargin = {x,y};
    if (nargin==4)
        varargin = [varargin,{Ttype}];
    end
end

% Calculate Affine-type transformation matrix:
inC = {'DefaultPixelValue',self.defVal};
if get(self.h.popup_T0type,'Value') == 6 % Warp
    np = min(size(self.points{1},1),size(self.points{2},1));
    off1 = self.cmiObj(1).img.voxsz.*self.cmiObj(1).img.dims(1:3)/2;
    off2 = self.cmiObj(2).img.voxsz.*self.cmiObj(2).img.dims(1:3)/2;
    p1 = ((self.points{1}(1:np,:)-repmat(off1([2,1,3]),np,1)))';
    p2 = ((self.points{2}(1:np,:)-repmat(off2([2,1,3]),np,1)))';
    M = reshape(p2([2,1,3],:)',1,[]);
    fM = reshape(p1([2,1,3],:)',1,[]);
    inC = [inC,{'FixedImageLandmarks',fM}];
else
    M = self.pts2M(varargin{:});
    M = [reshape(M(1:3,1:3),1,[]),M(1:3,4)'];
end

if ~isempty(M)
    % Set elxObj Initial Transform properties:
    self.elxObj.setTx0(M,...
            self.cmiObj(1).img.voxsz,...
            self.cmiObj(1).img.dims(1:3),...
            inC{:});
    set(self.h.checkbox_useExistingT,'Enable','on','Value',1);

    % Check for ouput directory:
    if isempty(self.odir)
        self.setOdir;
    end
    stat = exist(self.odir,'dir');

    % Call transformix:
    if stat
        tpfname = self.elxObj.saveTx0(fullfile(self.odir,...
            ['TransformParameters-snap-',datestr(now,'yyyymmdd'),'.txt']));
        fname = fullfile(self.odir,'elxtemp-in.mhd');
        stat = saveMHD(fname,self.cmiObj(2).img.mat(:,:,:,self.cmiObj(2).vec),...
                       [],self.cmiObj(2).img.dims(1:3).*self.cmiObj(2).img.voxsz);
        if stat
            tfxstr = self.elxObj.tfxCmd(self.odir,'tp',tpfname,'in',fname);

            % Run Transformix in new xterm window:
            namestr = 'Transformix';
            stat = ~system(['xterm -geometry 170x50 -T "',namestr,'"',...
                            ' -e ''',tfxstr{1},';''']);%csh''&']);

            % Append transformed image to Reference CMIobj
            if stat
                self.cmiObj(1).loadImg(true,fullfile(self.odir,'result.mhd'));
            end
        end
    end
end
    