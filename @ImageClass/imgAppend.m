% ImageClass function
function imgAppend(self,img,label)

if (nargin>1) && ~isempty(img)
    [d(1),d(2),d(3),d(4)] = size(img);
    if (nargin<3) || ~iscellstr(label) || (length(label)~=d(4))
        label = strcat({'Append-'},num2str((1:d(4))','%02u'))';
    end
    if self.check && all(self.dims(1:3)==d(1:3))
        self.check = true;
        ii = self.dims(4)+(1:d(4));
        self.mat(:,:,:,ii) = img;
        self.labels(ii) = label;
        self.scaleM(ii) = ones(1,d(4));
        self.scaleB(ii) = zeros(1,d(4));
        self.dims(4) = self.dims(4)+d(4);
        self.thresh(ii,:) = 1000000 * ones(d(4),1)*[-1,1];
        self.valExt(ii,:) = [squeeze(min(min(min(img)))),...
                             squeeze(max(max(max(img))))];
    elseif ~self.check
        self.mat = img;
        self.labels = labels;
        self.scaleM = ones(1,d(4));
        self.scaleB = zeros(1,d(4));
        self.dims = d;
        self.check = true;
        self.thresh = 1000000 * ones(d(4),1)*[-1,1];
        self.valExt = [ squeeze(min(min(min(img)))) ,...
                        squeeze(max(max(max(img)))) ];
    end
end