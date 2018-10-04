% CMIclass function

% Flip Image
function imgFlip(self,varargin)
if self.img.check
    if nargin==1 || isa(varargin{1},'matlab.ui.container.Menu')
        answer = inputdlg('Flip image in which dimension, X, Y, or Z (1, 2, 3)?');
        td = unique(str2num(char(answer)));
    elseif ismember(varargin{1},1:3)
        td = varargin{1};
    else
        error('Invalid input.');
    end
    self.img.flip(td);
    self.dispUDslice;
end