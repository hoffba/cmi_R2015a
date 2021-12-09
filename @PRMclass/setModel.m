% PRMclass function
function setModel(self,val)
% Set PRM options to default values for selected model

% Find saved PRM definitions
def_loc = fullfile(fileparts(which('cmi')),'PRMdefs');
prmdefs = dir(fullfile(def_loc,'PRMdef_*.mat'));
prmdefs = {prmdefs.name}';
prmlabels = extractBetween(prmdefs,'PRMdef_','.mat');

if (nargin == 1) || isempty(val)
    [val,ok] = listdlg('ListString',prmlabels,'SelectionMode','single',...
        'PromptString','Select a PRM model:');
    if ~ok, return; end
elseif ischar(val)
    val = find(strcmp(prmlabels,val),1);
end
if isempty(val)
    warning('Invalid input PRM model.');
else
    p = load(fullfile(def_loc,prmdefs{val}));
    self.setOpts('thresh',p.thresh,...
                 'cutoff',p.cutoff,...
                 'cmap',p.cmap,...
                 'prmmap',p.prmmap,...
                 'SPopts',p.SPopts);
end