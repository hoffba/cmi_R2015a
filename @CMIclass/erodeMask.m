% CMIclass function
% Erode mask
function erodeMask(self,~,~)
if self.img.mask.check
    prompt = {'Erosion Radius (pix):','Erosion Height (0:2D, >0:3D):'};
    def = {'3','0'}; % 3D dilation, x-y radius = 3, z radius = 1
    answer = str2double(inputdlg(prompt,'Erode VOI:',1,def));
    if ~isempty(answer) && ~any(isnan(answer))
        self.img.mask.morph('erode',answer);
        self.dispUDmask;
    end
end