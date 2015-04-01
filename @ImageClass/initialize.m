% ImageClass function
% re-initialize the object
function initialize(self)
self.mask.initialize;
self.prm.initialize;
self.mat = [];
self.dims = zeros(1,4);
self.voxsz = ones(1,3);
self.labels = {};
self.scaleM = [];
self.scaleB = [];
self.thresh = [];
self.prmBaseVec = 1;