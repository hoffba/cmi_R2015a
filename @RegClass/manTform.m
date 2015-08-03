% RegClass function
function M = manTform(self,x,~)
% Manually set initial transform

M = [];
tnames = {'Translation','Euler','Similarity','Affine'};
gchk = false; % GUI check
if ~self.cmiObj(1).img.check
    error('RegClass:manTform - Need to load a Reference Image before setting Tx0.');
elseif (nargin==1) || ((nargin==3) && ishandle(x))
    Ttype = get(self.h.popup_Transforms,'Value');
    gchk = true; x = [];
elseif isnumeric(x)
    if ismember(length(x),[3,6,7,12])
        self.elxObj.setTx0(x,self.cmiObj(1).img.voxsz,...
                           self.cmiObj(1).img.dims(1:3),...
                           'DefaultPixelValue',self.defVal);
    else
        error('RegClass:manTform - Invalid vector length for input parameters.')
    end
elseif ischar(x)
    Ttype = find(strcmp(x,tnames),1);
    x = [];
    if isempty(Ttype)
        error('RegClass:manTform - Invalid transform input.');
    else
        gchk = true;
    end
else
    error('RegClass:manTform - Invalid input.')
end

if gchk
    switch Ttype
        case 1
            wstr = 'No transform type has been selected.';
        case 2
            wstr = {'T_x','T_y','T_z'};
            def = zeros(1,3);
        case 3
            wstr = {'R_x','R_y','R_z','T_x','T_y','T_z'};
            def = zeros(1,6);
        case 4
            wstr = {'Q_x','Q_y','Q_z','T_x','T_y','T_z','S'};
            def = [zeros(1,6),1];
        case 5
            wstr = {'a_11','a_12','a_13','a_21','a_22','a_23',...
                    'a_31','a_32','a_33','T_x','T_y','T_z'};
            def = [1,0,0,0,1,0,0,0,1,0,0,0];
        case 6
            wstr = 'Manual inputs not available for warping transforms.';
        otherwise
            wstr = 'Invalid transform type.';
    end
    if ischar(wstr)
        warndlg(wstr);
    else
        % If transform exists, use values for default
        np = length(wstr);
        if ~isempty(self.elxObj.Tx0) && (length(self.elxObj.Tx0.TransformParameters)==np)
            def = self.elxObj.Tx0.TransformParameters;
        end
        def = cellfun(@num2str,num2cell(def),'UniformOutput',false);
        answer = inputdlg(wstr,['Manual Transform',tnames{Ttype-1}],1,def);
        x = cellfun(@str2double,answer);
        if any(isnan(x))
            error('Parameter inputs must be numeric.');
        end
        if Ttype==3
            x(1:3) = x(1:3)/180*pi; % degrees to radians
        end
    end
end

if ~isempty(x)
    M = self.par2affine(x(:)');
    self.elxObj.setTx0(reshape(M(1:3,:),1,[]),...
                       self.cmiObj(1).img.voxsz,...
                       self.cmiObj(1).img.dims(1:3),...
                       'DefaultPixelValue',self.T0defVal);
    self.showTx0(M(1:3,:));
    self.setTchk(true);
end
