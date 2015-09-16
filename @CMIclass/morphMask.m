% CMIclass function
% Perform morphological operation on mask (dilate/erode/open/close)
function morphMask(self,~,~)
if self.img.mask.check
    opts = {'d','dilate';'e','erode';'o','open';'c','close'};
    prompt = {'Operator ([d]ilate/[e]rode/[o]pen/[c]lose):',...
              'Dilation Radius (pix):',...
              'Dilation Height (0:2D, >0:3D):'};
    def = {'d','3','0'}; % 3D dilation, x-y radius = 3, z radius = 1
    answer = inputdlg(prompt,'Morph VOI:',1,def);
    if ~isempty(answer)
        o = max(strcmp(answer{1},opts),[],2);
        r = str2double(answer{2});
        h = str2double(answer{3});
        if ~(isnan(r)||isnan(h)) && any(o)
            self.img.mask.morph(opts{o,2},[r,h]);
            self.dispUDmask;
        end
    end
end