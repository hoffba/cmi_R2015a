% CMIclass function
% Set status of histogram
function setHistVis(self,x,~)
val = [];
if (nargin==3) && ishandle(x) && strcmp(get(x,'Tag'),'tools_disp_histo')
    val = ~strcmp(get(x,'Checked'),'on');
    if val
        str = 'on';
    else
        str = 'off';
    end
    set(x,'Checked',str);
elseif ~isempty(x) && islogical(x)
    val = x(1);
end
if islogical(val)
    self.histcheck = val;
    self.dispUDhist;
end