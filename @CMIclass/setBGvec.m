% CMIclass function
% Set background vector
function setBGvec(self,x,~)
val = [];
if (nargin==3) && ~isempty(x) && ishandle(x) && strcmp(get(x,'Tag'),'tools_disp_bgSel')
    val = str2double(inputdlg('Input image vector for background:',...
        'BG selection',1,{num2str(self.bgvec)}));
elseif ~isempty(x) && isnumeric(x)
    val = x(1);
end
if ~isempty(val) && ~isnan(val) && (val>0) && (val<=self.img.dims(4))
    self.bgvec = round(val);
    self.dispUDbg;
end