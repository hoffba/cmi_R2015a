% CMIclass function
% Set overlay check
function setOverCheck(self,x,~)
val = [];
if (nargin==3) && ishandle(x) && strcmp(get(x,'Tag'),'tools_disp_overlay')
    val = ~strcmp(get(x,'Checked'),'on');
elseif islogical(x) && ~isempty(x)
    val = x(1);
end
if ~isempty(val) && (self.overcheck ~= val)
    if val
        str = 'on';
    else
        str = 'off';
    end
    set(self.h.tools_disp_overlay,'Checked',str);
    self.overcheck = val;
    self.dispUDslice;
end