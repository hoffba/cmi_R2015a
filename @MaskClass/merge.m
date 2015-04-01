% MaskClass function
% Merge masks
function merge(self,meth,imask,slc,vdim)
% meth: 'replace' or 'add' or 'intersect' or 'remove'
% mask: 2D or 3D mask
% slc: (optional) specifies 2D slice
% vdim: (optional) specifies 2D slice dimension
opts = {'replace','add','intersect','remove','auto'};
if (nargin >= 3) && ~isempty(imask) && ismember(lower(meth),opts)
    tval = find(strcmpi(meth,opts));
    [d(1),d(2),d(3)] = size(imask);
    if ~self.check
        self.mat = false(self.dims(1:3));
    end
    if all(d==self.dims(1:3)) % Merge entire 3D mask
        slcheck = false;
        tmask = self.mat;
    elseif (nargin==5) && ismember(vdim,1:3) && ismember(slc,1:self.dims(vdim))
        slcheck = true;
        tmask = self.getSlice(vdim,1,slc);
    else
        slcheck = [];
    end
    
    % Merge the masks
    if ~isempty(slcheck)
        if isempty(tmask)
            tmask = zeros(d);
        end
        if tval == 5
            if self.check && any(tmask(:) & imask(:))
                tval = 4; % Remove
            else
                tval = 2; % Add
            end
        end
        switch tval
            case 1 %'replace'
                tmask = logical(imask);
            case 2 %'add'
                tmask = tmask | imask;
            case 3 %'intersect'
                tmask = tmask & imask;
            case 4 %'remove'
                tmask = tmask & ~imask;
        end
        if slcheck
            self.setSlice(tmask,vdim,1,slc)
        else
            self.mat = tmask;
        end
        self.check = any(self.mat(:));
    end
end