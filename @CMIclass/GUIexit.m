% CMIclass function
% exit the GUI and delete the CMIclass object
function GUIexit(self,~,~)
answer = questdlg('Would you really like to quit the program?','Quit CMI');
if strcmp(answer,'Yes')
    delete(self);
end


