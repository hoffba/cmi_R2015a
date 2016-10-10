% MaskClass function
% Load a new mask
function status = load(self,d,voxsz,fname,optsel)
% d = desired mask dimensions
% voxsz = voxel size of image containing mask
% fname = (cell) mask filename to load (incl. path)
status = false;
if (nargin >= 3) && ~isempty(d) && isnumeric(d) ...
        && (length(voxsz)==3) && isnumeric(voxsz)
    
    % If mask is already loaded, determine load/merge method
    ok = true;
    opts = {'replace','add','intersect','remove'};
    if self.check && ((nargin<5) || isempty(optsel))
        [sel,ok] = listdlg('ListString',opts,...
            'Name','VOI Merge Method','SelectionMode','single',...
            'PromptString','Merge masks?','ListSize',[160 100]);
    elseif isempty(optsel)
        sel = 1;
    else
        sel = find(strcmpi(optsel,opts),1);
        if isempty(sel)
            ok = false;
        end
    end
    
    if ok
        if ischar(fname)
            fname = {fname};
        end
        status = true;
        % Load new mask
        [tmask,~] = cmi_load(0,d,fname);
        % Perform mask math
        self.merge(opts{sel},tmask);
        end
    end
end