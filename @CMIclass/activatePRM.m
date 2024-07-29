% CMIclass function
% Activate PRM
function [chk,labels,vals] = activatePRM(self,x,~)
chk = [];
vals = [];
labels = {};

if ishghandle(x) && strcmp(x.Tag,'analysis_PRM')
    chk = ~strcmp(get(self.h.analysis_PRM,'Checked'),'on');
elseif ~isempty(x) && islogical(x)
    chk = x(1);
end

if ~isempty(chk)
    if (chk ~= self.prmcheck)
        self.prmcheck = chk;
        if chk && self.img.mask.check
            [labels,vals] = self.img.calcPRM(self.vec);
        end
        if self.guicheck
            str = {'Off','On'};
            if chk
                str = str([2,1]);
            end
            set(self.h.analysis_PRM,'Checked',lower(str{1}));
            h = findall(self.hfig,'ToolTipString','Data Cursor');
            set(h,'Visible',str{2});
            self.dispUDslice;
            self.dispUDcmap;
        end
    else
        chk = self.prmcheck;
    end
end