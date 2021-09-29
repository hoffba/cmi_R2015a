% ElxClass function
function stat = setFormat(self,str)
% Set outfmt parameter: format for saving images

stat = false;
if ischar(str) && ismember(str,{'.mhd','.nii','.nii.gz'})
    self.outfmt = str;
end
