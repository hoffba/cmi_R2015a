% ImageClass function
% Set image values inside/outside mask to specific value
function setMaskVal(self,vec,mval,ival)
if (nargin==4) && isnumeric(vec) && (isnumeric(mval) || islogical(mval)) ...
        && isnumeric(ival) && ~(isempty(vec) || isempty(mval) || isempty(ival))
    nv = length(vec);
    if length(ival)~=nv
        ival = ival(1)*ones(1,nv);
    end
    vec = round(vec);
    mval = logical(mval(1));
    ival = double(ival);
    d4 = self.dims(4);
    if all((vec>0) & (vec<=d4)) && ~any(isnan(ival))
        for i = 1:nv
            timg = self.mat(:,:,:,i);
            timg(self.mask.mat==mval) = ival(i);
            self.mat(:,:,:,i) = timg;
        end
    else
        error(['Invalid image index: ',num2str(vec)]);
    end
else
    error('Invalid inputs.');
end