function addElxStep(self, type_str, varargin)
% Add step to the Elastix Schedule

opts = {'Translation','Euler','Similarity','Affine','Warp'};

if nargin<1
    [sel,ok] = listdlg('ListString',opts,...
        'SelectionMode','single',...
        'Name','Transform Type');
    if ok
        type_str = opts{sel};
    else
        type_str = '';
    end
end

if ischar(type_str) && ~isempty(type_str)
    if ismember(type_str,opts)
        stat = self.elxObj.addStep(type_str);
        if stat
            % Set additional parameters:
            ind = length(self.elxObj.Schedule);     
            if nargin>1
                self.setElxPar(ind,varargin{:});
            end
            
            % Update GUI objects:
            if self.guicheck
                if strcmp(get(self.h.listbox_Tforms,'Enable'),'off')
                    set(self.h.listbox_Tforms,'Enable','on');
                end
                str = [get(self.h.listbox_Tforms,'String');{type_str}];
                set(self.h.listbox_Tforms,'String',str);
                self.selectTform(length(str));
            end
        end
    else
        error('Invalid transform type: %s',type_str);
    end
end