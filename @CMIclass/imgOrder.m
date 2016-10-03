% CMIclass method
function imgOrder(self,x,~)
% Reorders images (4D)

if self.img.check
    if nargin==1 || isa(x,'matlab.ui.container.Menu') % GUI for input:
        x = inputdlg(sprintf('New image order (1-%u):',self.img.dims(4)),...
            'Image Order',1,{num2str(1:self.img.dims(4))});
        if isempty(x)
            return;
        else
            x = str2num(x{1});
        end
    end
    stat = self.img.reorder(x);
    if stat
        self.dispUDslice;
    end
end