% RegClass function
% Switch for the radio buttons and point selection enabeling
function radSwitch(self, hObject, ~)
if ishandle(hObject)
    h = [ self.h.radRef  , self.h.radHom ];
    self.imgSelect = (hObject==h(2))+1;
    if self.imgSelect==2
        h = h([2,1]);
    end
    set(h(1),'Value',1);
    set(h(2),'Value',0);
elseif ismember(hObject,[1,2])
    self.imgSelect = hObject;
end