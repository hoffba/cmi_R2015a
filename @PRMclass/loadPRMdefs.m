% PRMclass function
% Load PRM definitions from .mat file
function stat = loadPRMdefs(self,fname)

stat = false;
if nargin==1
    [fname,fdir] = uigetfile('*PRMdef*.mat','Load PRM Defs');
    if ischar(fname)
        fname = fullfile(fdir,fname);
    end
end

if ischar(fname) && exist(fname,'file')
    t = load(fname);
    t = [fieldnames(t),struct2cell(t)]';
    self.setOpts(t{:});
    stat = true;
end