function str = sysCmd(self,odir,varargin)

str = {};
mxCores = feature('numCores');

p = inputParser;
p.addParameter('title','Elastix / Transformix',@ischar);
p.addRequired('odir',@(x)ischar(x)&&exist(x,'dir'));
p.addParameter('wait',false,@islogical);
p.addParameter('cleanup',false,@islogical);
% Elastix
p.addParameter('f','',@ischar);
p.addParameter('m','',@ischar);
p.addParameter('fMask','',@ischar);
p.addParameter('mMask','',@ischar);
p.addParameter('p','',@(x)ischar(x)||iscellstr(x));
    % (p can be explicitly set or saved from GUI inputs)
p.addParameter('t0','',@ischar);
p.addParameter('threads',[],@(x)isnumeric(x)&&(x>0)&&(x<=mxCores));
% Transformix:
p.addParameter('tp','',@(x)ischar(x)||iscellstr(x));
p.addParameter('in','',@(x)ischar(x)||iscellstr(x));
p.addParameter('jac',false,@islogical);
p.addParameter('jacmat',false,@islogical);
p.addParameter('def',false,@islogical);
p.addParameter('outfn',{},@(x)ischar(x)||iscellstr(x));
% Parse inputs:
parse(p,odir,varargin{:});
iflds = setdiff(p.Parameters,p.UsingDefaults);
p = p.Results;

% Need to remove trailing "\" on directories in Windows:
if strcmp(p.odir(end),filesep)
    p.odir(end) = [];
end

% Generate Elastix call:
if all(ismember({'f','m'},iflds))
    str = {self.elxCmd(p.f,...
                       p.m,...
                       p.odir,...
               'fMask',p.fMask,...
               'mMask',p.mMask,...
             'threads',p.threads)};
end
% Generate Transformix call:
ii = ismember({'tp','in','jac','jacmat','def'},iflds);
if ii(1) && any(ii(2:end))
    str = [str;self.tfxCmd(p.odir,...
                        'in',p.in,...
                     'outfn',p.outfn,...
                        'tp',p.tp,...
                       'jac',p.jac,...
                    'jacmat',p.jacmat,...
                       'def',p.def)];
end
if ~isempty(str)
    % Cleanup string:
    if p.cleanup
        if ispc
            tstr = {}; % not programmed yet
        else
            tstr = {['find ',p.odir,' -name "elxtemp-*" -exec rm -f {} \;']};
        end
        str = [str;tstr];
    end
    % Option to keep terminal window open:
    if ~p.wait
        if ispc
            tstr = {};
        else
            tstr = 'csh';
        end
        str = [str;tstr];
        wstr = '&';
    else
        wstr = '';
    end
    % Concatenate system calls:
    str = sprintf(['%s',repmat([self.sepstr,'%s'],1,numel(str)-1)],str{:});
    if ispc
        str = [self.xtstr,'title ',p.title,self.sepstr,str,wstr];
    else
        str = [self.xtstr,'-T "',p.title,'" -e ''',str,'''',wstr];
    end
end
