% CMIcalss function
% Activate PRM
function val = activatePRM(self,x,~)
val = [];
h = self.h.analysis_PRM;

if (nargin==3) && (x==h) %ishandle(x) && strcmp(get(x,'Tag'),'analysis_PRM')
    val = ~strcmp(get(h,'Checked'),'on');
elseif ~isempty(x) && islogical(x)
    val = x(1);
end

if ~isempty(val) && islogical(val)
    % Set GUI tool display check
    if val
        str = 'on';
    else
        str = 'off';
    end
    set(h,'Checked',str);
    
    % Figure tools opposite of check
    if val
        str = 'Off';
    else
        str = 'On';
    end
    
    if (val ~= self.prmcheck)
        self.prmcheck = val;
        if val && self.img.mask.check
            self.img.calcPRM(self.vec);
        end
        h = findall(self.hfig,'ToolTipString','Data Cursor');
        set(h,'Visible',str);
        self.dispUDslice;
        self.dispUDcmap;
    else
        val = self.prmcheck;
    end
end