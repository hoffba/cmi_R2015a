% ImageClass function
% Load image
function status = loadImg(self,appendcheck,fname)
if (nargin <2)
    appendcheck = false;
end
if (nargin < 3)
    fname = {};
end
if ischar(fname)
    fname = {fname};
end
d = [];
if appendcheck
    d = self.dims(1:3);
end
[img,label,fov,fnameOut] = cmi_load(1,d,fname);
% img(isnan(img)) = 0;
if ~isempty(img)
%     if size(img,3)==1
%         img = permute(img,[1 2 4 3]);
%         label = label(1);
%         fov(3) = mean(fov(1:2));
%     end

    dd = size(img);
    if length(dd)<3
        dd(3) = 1;
    end
    
    % If reading in a complex image, separate into real and imaginary
    % components
    if ~isreal(img)
        nv = size(img,4);
        img = cat(4,real(img),imag(img));
        label = strcat([repmat({'Re-'},[1,nv]),repmat({'Im-'},[1,nv])],...
                    label);
    end
    
    if appendcheck && self.check && all(self.dims(1:3)==dd(1:3))
        nnv = size(img,4);
        self.mat = cat(4,self.mat,img);
        self.dims(4) = size(self.mat,4);
        self.labels = [self.labels label];
        self.thresh = [self.thresh; ([-1 1]'*inf(1,nnv))'];
        self.scaleM = [self.scaleM ones(1,nnv)];
        self.scaleB = [self.scaleB zeros(1,nnv)];
        self.valExt = [self.valExt ; [squeeze(min(min(min(img,[],1),[],2),[],3)),...
                                      squeeze(max(max(max(img,[],1),[],2),[],3))]];
    else
        self.prmBaseVec = 1;
        [self.dir,self.name] = fileparts(fnameOut);
        [d(1),d(2),d(3),d(4)] = size(img);
        self.mat = img;
        self.labels = label;
        self.voxsz = fov./d(1:3);
        self.dims = d;
        self.mask.setDims(d(1:3));
        self.thresh = ([-1 1]'*inf(1,d(4)))';
        self.scaleM = ones(1,d(4));
        self.scaleB = zeros(1,d(4));
        self.valExt = [squeeze(min(min(min(img,[],1),[],2),[],3)),...
                       squeeze(max(max(max(img,[],1),[],2),[],3))];
        self.mask.initialize(d(1:3));
        self.prm.initialize;
    end
    self.check = true; % Image is now available
    status = true; % Load completed successfully
else
    status = false; % Nothing was loaded
end