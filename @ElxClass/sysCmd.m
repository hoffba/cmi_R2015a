function str = sysCmd(self,odir,varargin)

str = '';
istr = {};
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
p.addParameter('tp','',@ischar);
p.addParameter('in','',@(x)ischar(x)||iscellstr(x));
p.addParameter('jac',false,@islogical);
p.addParameter('jacmat',false,@islogical);
p.addParameter('def',false,@islogical);
p.addParameter('outfn','',@ischar);
% Parse inputs:
parse(p,odir,varargin{:});
iflds = setdiff(p.Parameters,p.UsingDefaults);

% Generate Elastix call:
if all(ismember({'f','m'},iflds))
    istr = self.elxCmd(p.Results.f,...
                       p.Results.m,...
                       p.Results.odir,...
               'fMask',p.Results.fMask,...
               'mMask',p.Results.mMask,...
             'threads',p.Results.threads);
end
% Generate Transformix call:
ii = ismember({'tp','in','jac','jacmat','def'},iflds);
if ii(1) && any(ii(2:end))
    istr = [istr,self.tfxCmd(p.Results.odir,...
                        'in',p.Results.in,...
                     'outfn',p.Results.outfn,...
                        'tp',p.Results.tp,...
                       'jac',p.Results.jac,...
                    'jacmat',p.Results.jacmat,...
                       'def',p.Results.def)];
end

if ~isempty(istr)
    % Generate executable:
    str = istr{1};
    if length(istr)>1
        istr = sprintf([self.sepstr,' %s'],istr{2:end});
    else
        istr = '';
    end
    istr = [str,istr];
    % Cleanup string:
    if p.Results.cleanup
        if ispc
            custr = '';
        else
            custr = ['find ',p.Results.odir,' -name "elxtemp-*" -exec rm -f {} \;'];
        end
        istr = [istr,self.sepstr,custr];
    end
    % Generate system call:
    str = self.xtstr;
    if ispc
        str = [str,'title ',p.Results.title,self.sepstr,istr];
    else
        str = [str,'-T "',p.Results.title,'" -e ''',istr];
    end
    % Option to wait for execution:
    if ~p.Results.wait
        if ispc
            waitstr = '&';
        else
            waitstr = ' ;csh;''&';
        end
    elseif ~ispc
        waitstr = '''';
    end
    str = [str,waitstr];
end

