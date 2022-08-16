function result = proc_unreg(self,ei_tag)
% Process unregistered analysis

result = [];
if ismember('exp',ei_tag)
    [img,label,ATmap] = self.getData('ct_exp','seg_exp','scatnet');
    if ~(isempty(img) || isempty(label))
        result = CTlung_Unreg('exp',img,self.dat.ct_exp.info.PixelDimensions,label,ATmap);
    end
end
if ismember('ins',ei_tag)
    [img,label] = self.getData('ct_exp','seg_exp');
    if ~(isempty(img) || isempty(label))
        res = CTlung_Unreg('ins',img,self.dat.ct_ins.info.PixelDimensions,label);
        result = mergeStructs(result,res);
    end
end