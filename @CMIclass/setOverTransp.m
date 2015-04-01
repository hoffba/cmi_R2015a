% CMIclass function
function setOverTransp(self,~,~)
% Set overlay transparency
answer = cellfun(@str2double,...
                 inputdlg('Input image transparency (0:1):',...
                          'Transparency',1,{num2str(self.dalpha)}));
if ~isempty(answer) && (answer>0) && (answer<=1)
    self.dalpha = answer;
    self.dispUDimg;
end