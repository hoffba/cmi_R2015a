% CMIcalss function
% Generates isosurface on image or mask 3D dataset
function isosurf(self,~,~)
if self.img.check
    ok = true;
    tvec = self.vec;
    if self.img.mask.check
        answer = questdlg('Show isosurface for ...','Isosurface','Image','Mask','Cancel','Image');
        if strcmp(answer,'Mask')
            tvec = 0;
            val = 0.5;
        elseif strcmp(answer,'Cancel')
            ok = false;
        end
    end
    if tvec
        val = str2double(inputdlg('Input threshold for surface:','Isosurface'));
    end
    if ok && ~isempty(val) && ~isnan(val)
        self.img.isosurf(tvec,val);
    end
end