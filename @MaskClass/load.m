% MaskClass function
% Load a new mask
function status = load(self,d,voxsz,fname)
% d = desired mask dimensions
% voxsz = voxel size of image containing mask
% fname = (cell) mask filename to load (incl. path)
status = false;
if (nargin >= 3) && ~isempty(d) && isnumeric(d) ...
        && (length(voxsz)==3) && isnumeric(voxsz)
    % If mask is already loaded, determine load/merge method
    opts = {'replace','add','intersect','remove'};
    if self.check && isempty(fname) && ~any(d ~= self.dims(1:3))
        [sel,ok] = listdlg('ListString',opts,...
            'Name','VOI Merge Method','SelectionMode','single',...
            'PromptString','Merge masks?','ListSize',[160 100]);
    else
        ok = 1;
        sel = 1;
    end
    if ok
        if ischar(fname)
            fname = {fname};
        end
        status = true;
        % Load new mask
        [tmask,~] = cmi_load(0,d,fname);
        %tmask = permute(tmask,self.calcDimOrder);
%         for tf = 1:3
%             if self.transp(tf)
%                 tmask = flipdim(img,tf);
%             end
%         end
        % Perform mask math
        self.merge(opts{sel},tmask);
    end
end