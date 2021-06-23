% CMIclass function
% Activate PRM
function val = activatePRM(self,x,~)
val = [];

if ishghandle(x) && strcmp(x.Tag,'analysis_PRM')
    val = ~strcmp(get(self.h.analysis_PRM,'Checked'),'on');
elseif ~isempty(x) && islogical(x)
    val = x(1);
end

if ~isempty(val)
    if (val ~= self.prmcheck)
        self.prmcheck = val;
        if val && self.img.mask.check
            self.img.calcPRM(self.vec);
        end
        if self.guicheck
            str = {'Off','On'};
            if val
                str = str([2,1]);
            end
            set(self.h.analysis_PRM,'Checked',lower(str{1}));
            h = findall(self.hfig,'ToolTipString','Data Cursor');
            set(h,'Visible',str{2});
            self.dispUDslice;
            self.dispUDcmap;
        end
    else
        val = self.prmcheck;
    end
end