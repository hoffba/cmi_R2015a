% PRMclass function
function stat = setOpts(self,varargin)
% Set PRM options -- Name/Value pairs
%           varargin = {'Name',Value,...} --> manual setting

stat = false;
p = inputParser;
p.KeepUnmatched = true;
p.addOptional('thresh',nan,@(x) isnumeric(x) && (size(x,2)==4));
p.addOptional('cutoff',nan,@(x) isnumeric(x) && any(size(x,2)==[0,3]));
p.addOptional('prmmap',nan,@(x) iscell(x) && (size(x,2)==2) ...
                && all(cellfun(@islogical,x(:,1))) ...
                && iscellstr(x(:,2)));
p.addOptional('cmap'  ,nan,@(x) isnumeric(x) && (size(x,2)==3));
p.addOptional('labels',nan,@(x) iscellstr(x));
p.addOptional('filtchk',nan,@(x) islogical(x) && length(x)==1);
p.addOptional('filttype',nan,@(x) ischar(x));
p.addOptional('filtstr',nan,@(x) isvector(x) && length(x)==2);
p.addOptional('SPopts',nan,@isstruct);
p.addOptional('guicheck',nan,@isscalar);
p.parse(varargin{:});
opts = rmfield(p.Results,p.UsingDefaults);

if isfield(opts,'thresh')
    if isempty(opts.thresh)
        tthresh = zeros(0,size(self.thresh,2));
    else
        tthresh = opts.thresh;
    end
    tthresh(:,1:2) = round(tthresh(:,1:2)); % dimensions must be integer
    ind = (tthresh(:,1)~=tthresh(:,2)) & ~any(tthresh(:,1:2)<0,2);
    tthresh = tthresh(ind',:); % only use valid thresholds
    self.thresh = tthresh;
    % Extract relevent image indices:
    mvals = mode(self.thresh(:,1:2)); % for the scatterplot
    tval = unique(self.thresh(:,1:2));
    tval(ismember(tval,mvals)) = [];
    self.dvec = [mvals,tval'];
    
    stat = true;
end

if isfield(opts,'cutoff')
    tcutoff = opts.cutoff;
    if isempty(tcutoff)
        self.cutoff = zeros(0,3);
    else
        ind = (tcutoff(:,1)>=0);
        tcutoff = tcutoff(ind',:);
        if isempty(tcutoff)
            tcutoff = zeros(0,3);
        else
            tcutoff(:,2:3) = sort(tcutoff(:,2:3),2); % sort for [min max] order
        end
        self.cutoff = tcutoff;
    end
    stat = true;
end 

% prmmap logical checks must have ncol == nthresh
nth = size(self.thresh,1);
if isfield(opts,'prmmap') && ...
        all(cellfun(@(x)size(x,2),opts.prmmap(:,1)) == nth)
    self.prmmap = opts.prmmap;
    self.nprm = size(self.prmmap,1);
    stat = true;
end

if isfield(opts,'cmap') && (size(p.Results.cmap,1)==size(self.prmmap,1))
    self.cmap = opts.cmap;
    stat = true;
end

if isfield(opts,'labels') && (length(p.Results.labels)==length(self.dvec))
    self.dlabels = opts.labels;
    stat = true;
end

if isfield(opts,'filtchk')
    self.filtchk = opts.filtchk;
    stat = true;
end

ftypes = {'wiener','gauss','median','average'};
if isfield(opts,'filttype') && ismember(lower(opts.filttype),ftypes)
    self.filttype = opts.filttype;
    stat = true;
end

if isfield(opts,'filtstr')
    self.filtstr = opts.filtstr;
    stat = true;
end

if isfield(opts,'SPopts')
    for fldn = fieldnames(self.SPopts)'
        if isfield(opts.SPopts,fldn{1})
            self.SPopts.(fldn{1}) = opts.SPopts.(fldn{1});
        end
    end
end
        
if isfield(opts,'guicheck')
    self.guicheck = opts.guicheck;
    stat = true;
end
