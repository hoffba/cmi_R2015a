% CMIclass function
% Toggle checks for analysis model selection
function modCheck(self,ht,~)
if (nargin==3) && ishandle(ht)
    thandles = [self.h.analysis_genmod,...
                self.h.analysis_perfmod,...
                self.h.analysis_diffmod     ];
    ind = (ht==thandles);
    set(thandles(ind),'Checked','on');
    set(thandles(~ind),'Checked','off');
    self.img.setModelType(find(ind,1));
    set(self.h.analysis_modtype,'Label',...
        ['Model Type: ',self.img.model.typeStr])
    set(self.h.analysis_modsel,'Label',...
        ['Model: ',self.img.model.getModName]);
    if isa(self.img.model,'DCEclass')
        set([self.h.analysis_dceopts,self.h.analysis_setAIF],'Enable','on')
    else
        set([self.h.analysis_dceopts,self.h.analysis_setAIF],'Enable','off')
    end
end