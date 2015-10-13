% CMIclass function
% Allows user to manually translate images
function manreg(self)
if self.img.check
    
    v = 1;
    if self.img.dims(4)>1
        answer = questdlg('Move all images or only the current?',...
            'Shift Image','All','Current','Current');
        
    end
    
end