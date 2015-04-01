% RegClass function
function delPoints(self,hObject,~)
% Delete selected points from list

iimg = find(strcmp(hObject.Tag(end-2:end),{'Ref','Hom'}),1);
if ~isempty(self.points{iimg})
    if iimg == 1
        h = self.h.list_ref;
    else
        h = self.h.list_hom;
    end
    p = self.points{iimg};
    i = get(h,'Value');
    if ~isempty(i)
        p(i,:) = [];
        self.points{iimg} = p;
        self.plotPts(iimg);
    end
end