% MaskClass function
% Save current mask
function status = save(self,voxsz,fname)
if self.check
    if (nargin == 1) || ~isnumeric(voxsz) || (length(voxsz)~=3)
        voxsz = ones(1,3);
    end
    if (nargin == 3) && ischar(fname)
        [path,name,ext] = fileparts(fname);
        if exist(path,'dir')
            if isempty(ext) && ~isempty(name) && exist(path,'dir')
                % Default if no extension was input
                fname = fullfile(path,name,'.mask'); 
            elseif any(strcmpi(ext,{'.mask','.fld'}))
                % Mask only saved as MASK or FLD file type
                fname = fullfile(path,name,ext); 
            end
        end
    else
        fname = '';
    end
%     tmask = self.mat(:,:,:);
    fov = self.dims(1:3).*voxsz;
    if (nargin<2) % fname was not input
        fname = [];
    end
    info = struct('dircos',self.dircos,'slcpos',self.slcpos);
    status = cmi_save(logical(self.mat),{'VOI'},fov,fname,info);
else
    status = false;
end