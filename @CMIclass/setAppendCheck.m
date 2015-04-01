% CMIclass function
% Set append check
function setAppendCheck(self,x,~)
val = [];
if nargin>1
    if (nargin==3) && ishandle(x) && strcmp(get(x,'Tag'),'file_append')
        val = ~strcmp(get(x,'Checked'),'on');
        if val % if off, turn on
            str = 'on';
        else % if on, turn off
            str = 'off';
        end
        set(x,'Checked',str);
    elseif ~isempty(x) && isnumeric(x)
        val = logical(x(1));
    end
end
if ~isempty(val)
    self.appendcheck = val;
end