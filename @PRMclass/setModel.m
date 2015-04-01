% PRMclass function
function setModel(self,val)
% Set PRM options to default values for selected model
if (nargin == 1) || isempty(val) || ~isnumeric(val)...
        || (val(1)>length(self.prmdef)) || (val(1)<1)
    opts = squeeze(struct2cell(self.prmdef));
    [val,ok] = listdlg('ListString',opts(1,:),'SelectionMode','single',...
        'PromptString','Select a PRM model:');
else
    val = round(val);
    ok = true;
end
if ok
    self.thresh = self.prmdef(val).thresh;
    self.cutoff = self.prmdef(val).cutoff;
    self.cmap = self.prmdef(val).cmap;
    self.prmmap = self.prmdef(val).prmmap;
    self.nprm = size(self.prmmap,1);
    self.dvec = unique(self.thresh(:,1:2));
    self.dlabels = {'Expiration','Inspiration'};
end