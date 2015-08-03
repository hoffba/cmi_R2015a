% ImageClass function
% Set image values inside/outside mask to specific value
function setMaskVal(self,vec,mval,ival)
if (nargin==4) && isnumeric(vec) && (isnumeric(mval) || islogical(mval)) ...
        && isnumeric(ival) && ~(isempty(vec) || isempty(mval) || isempty(ival))
    vec = round(vec);
    mval = logical(mval(1));
    ival = double(ival(1));
    d4 = self.dims(4);
    if all((vec>0) & (vec<=d4)) && ~(isnan(vec) || isnan(ival))
        tmask = repmat(self.mask.mat,1,1,1,d4);
        tmask(:,:,:,~ismember(1:d4,vec)) = false;
        self.mat(tmask==mval) = ival;
    else
        error(['Invalid image index: ',num2str(vec)]);
    end
else
    error('Invalid inputs.');
end