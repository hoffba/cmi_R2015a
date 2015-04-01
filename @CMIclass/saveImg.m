% CMIclass function
% Save image
function saveImg(self,x,~)
if self.img.check
    fname = '';
    if ischar(x)
        fname = x;
    end
    fname = self.img.saveImg([],fname,self.vec);
    if fname
        set(self.hfig,'Name',self.img.name)
    end
end