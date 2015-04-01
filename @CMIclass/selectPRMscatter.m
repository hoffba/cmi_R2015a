% CMIclass function
function selectPRMscatter(self,~,~)

str = get(self.h.analysis_prmselect,'Checked');
if strcmp(str,'on')
    str = 'off';
    stat = true;
else
    str = 'on';
    stat = self.img.selectPRMscatter();
end
if stat % Update PRM overlay
    set(self.h.analysis_prmselect,'Checked',str);
    self.dispUDslice;
end
