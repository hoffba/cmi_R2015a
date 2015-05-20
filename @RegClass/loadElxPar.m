% RegClass function
function loadElxPar(self,fname,~)
% loadElxPar
% loadElxPar('path/fname.txt')
% loadElxPar({'fname1','fname2',...})

if (nargin==1) || ~ischar(fname)
    fname = {};
end
stat = self.elxObj.loadPar(fname);
if stat
    ns = length(self.elxObj.Schedule);
    % Update list of Elastix steps:
    str = cell(1,ns);
    for i = 1:ns
        tstr = self.elxObj.getPar(i,'Transform');
        ind = strfind(tstr,'Transform');
        tstr = tstr(1:ind-1);
        if strcmp(tstr,'BSpline')
            tstr = 'Warp';
        end
        str{i} = tstr;
    end
    set(self.h.listbox_Tforms,'String',str,'Enable','on');
    % Snap to last step on the list and update GUI:
    self.selectTform(ns);
end
