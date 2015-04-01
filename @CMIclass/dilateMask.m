% CMIclass function
% Dilate mask
function dilateMask(self,~,~)
if self.img.mask.check
    prompt = {'Dilation Radius (pix):','Dilation Height (0:2D, >0:3D):'};
    def = {'3','0'}; % 3D dilation, x-y radius = 3, z radius = 1
    answer = str2double(inputdlg(prompt,'Dilate VOI:',1,def));
    if ~isempty(answer) && ~any(isnan(answer))
        self.img.mask.morph('dilate',answer);
        self.dispUDmask;
    end
end