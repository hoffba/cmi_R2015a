% CMIclass function
% Set CMIclass property
function setPlotProp(self,varargin)
% Set CMIclass property value

propns = {'voispec','voimarksz','thspec','thmarksz'};

cspecs = {'y','m','c','r','g','b','w','k'};
mspecs = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};

if (nargin==3) && isa(varargin{1},'matlab.ui.control.UIControl')
    h = varargin{1};
    N = propns(strcmp(h.Tag,{'edit_voiLSpec','edit_voiMkSz',...
                             'edit_thrLSpec','edit_thrMkSz'}));
    V = {h.String};
elseif length(varargin)>1
    N = varargin(1:2:end);
    V = varargin(2:2:end);
else
    N = {};
end

nv = length(N);
if (nv>0) && (length(V)==nv)
    for i = 1:nv
        switch N{i}
            case 'voispec'
                if ischar(V{i})
                    n = length(V{i});
                    if (n>0) && ismember(V{i}(1),mspecs)
                        self.voispec(1) = V{i}(1);
                    end
                    if (n>1) && ismember(V{i}(2),cspecs)
                        self.voispec(2) = V{i}(2);
                    end
                    set(self.h.edit_voiLSpec,'String',self.voispec);
                    self.dispUDroi;
                end
            case 'voimarksz'
                if ischar(V{i})
                    V{i} = str2double(V{i});
                end
                if isnumeric(V{i}) && ~isnan(V{i})
                    self.voimarksz = V{i};
                end
                set(self.h.edit_voiMkSz,'String',num2str(self.voimarksz));
                self.dispUDroi;
            case 'thspec'
                if ischar(V{i})
                    n = length(V{i});
                    if (n>0) && ismember(V{i}(1),mspecs)
                        self.thspec(1) = V{i}(1);
                    end
                    if (n>1) && ismember(V{i}(2),cspecs)
                        self.thspec(2) = V{i}(2);
                    end
                    set(self.h.edit_thrLSpec,'String',self.thspec);
                    self.dispUDthreshplot;
                end
            case 'thmarksz'
                if ischar(V{i})
                    V{i} = str2double(V{i});
                end
                if isnumeric(V{i}) && ~isnan(V{i})
                    self.thmarksz = V{i};
                end
                set(self.h.edit_thrMkSz,'String',num2str(self.thmarksz));
                self.dispUDthreshplot;
        end
    end
    self.dispUDroi;
end