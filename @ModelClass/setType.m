% ModelClass function
function setType(self,ind)
if (nargin==2)
    if any(ind==(1:length(self.typestrs)))
        self.mtype = ind;
    elseif ischar(ind)
        tind = find(strcmp(ind,self.typestrs),1);
        if ~isempty(tind)
            self.mtype = tind;
        end
    end
end