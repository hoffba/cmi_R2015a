% CMIclass function

% Flip Image
function imgFlip(self,~,~)
if self.img.check
    answer = inputdlg('Flip image in which dimension, Y, X, or Z (1, 2, 3)?');
    td = unique(str2num(char(answer)));
    if ~isempty(td)
        self.img.flip(td);
        self.dispUDslice;
    end
end