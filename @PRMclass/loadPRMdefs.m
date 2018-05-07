% PRMclass function
% Load PRM definitions from .mat file
function stat = loadPRMdefs(self,fname)

stat = false;
if nargin==1
    [fname,fdir] = uigetfile('*PRMdef*.mat','Load PRM Defs');
    if ~ischar(fname)
        return;
    end
else
    [fdir,fname]  = fileparts(fname);
    if isempty(fdir)
        fdir = fullfile(fileparts(which('cmi')),'PRMdefs');
    end
    fname = [fname,'.mat'];
end
fname = fullfile(fdir,fname);
if ischar(fname) && exist(fname,'file')
    t = load(fname);
    t = [fieldnames(t),struct2cell(t)]';
    self.setOpts(t{:});
    stat = true;
end