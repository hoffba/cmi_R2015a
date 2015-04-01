% CMIclass function
% Set PRM options
function setPRMopts(self,~,~)
stat = self.img.setPRMopts(self.vec); % no inputs -> prompt user for input
if stat
    if self.img.prm.check && self.prmcheck
        self.img.calcPRM(self.vec);
        self.dispUDcmap;
    end
    self.dispUDimg;
end