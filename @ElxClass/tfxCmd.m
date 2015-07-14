% ElxClass function
% Generate system call to elastix
function str = tfxCmd(self,odir,varargin)
% Inputs: odir  = output directory
%         Name/Value pairs (need at least one of the following:)
%           'in'    / Images to transform (.mhd)
%           'outfn' / File names for transformed 'in' images
%           'jac'   / Flag to generate Jacobian map (TF)
%           'def'   / Flag to generate deformation maps (TF)
%           'jacmat'/ Flag to generate Jacobian matrices (TF)
%         Name/Value pairs (optional)
%           'tp'    / TransformParameterFile (.txt)
%               * if absent, uses ElxClass.Tx0 property

str = '';
% Parse variable inputs:
p = inputParser;
addRequired(p,'odir',@isdir);
addParameter(p,'in','',@(x)ischar(x)||iscellstr(x));
addParameter(p,'outfn',{},@(x)ischar(x)||iscellstr(x));
addParameter(p,'tp','',@(x)ischar(x)||iscellstr(x));
addParameter(p,'jac',false,@islogical);
addParameter(p,'jacmat',false,@islogical);
addParameter(p,'def',false,@islogical);
parse(p,odir,varargin{:})
pp = p.Results;

% Determine tranform parameters:
if isempty(pp.tp) && ~isempty(self.Tx0)
    % If no file was input, save current Tx0 in ExlClass object
    pp.tp = self.saveTx0(fullfile(pp.odir,...
                ['TransformParameters-',datestr(now,'yyyymmdd'),'.txt']));
end
if ischar(pp.in)
    pp.in = {pp.in};
end
if ischar(pp.outfn)
    pp.outfn = {pp.outfn};
end
if ischar(pp.tp)
    pp.tp = {pp.tp};
end
nf = length(pp.in);
if length(pp.tp)~=nf
    pp.tp = repmat(pp.tp,nf,1);
end

if ~isempty(pp.tp) && (~isempty(pp.in) || pp.jac || pp.jacmat || pp.def)
    outchk = length(pp.outfn)==nf;
    str = cell(1+2*outchk,nf);
    for i = 1:nf
        str{1,i} = [fullfile(self.elxdir,'transformix'),...
               ' -out "',pp.odir,'"',...
               ' -tp "',pp.tp{i},'"'];
        if ~isempty(pp.in)
            str{1,i} = [str{1,i},' -in "',pp.in{i},'"'];
        end
        if pp.jac
            str{1,i} = [str{1,i},' -jac all'];
            pp.jac = false; % only need to perform once
        end
        if pp.def
            str{1,i} = [str{1,i},' -def all'];
            pp.def = false; % only need to perform once
        end
        if pp.jacmat
            str{1,i} = [str{1,i},' -jacmat all'];
            pp.jacmat = false; % only need to perform once
        end
        if outchk && ~isempty(pp.outfn) && ~isempty(pp.outfn{i})
            if ispc
                cpexec = 'copy';
            else
                cpexec = 'cp';
            end
            str{2,i} = [cpexec,' "',fullfile(pp.odir,'result.mhd'),'" "',...
                        pp.outfn{i},'"'];
            str{3,i} = [cpexec,' "',fullfile(pp.odir,'result.raw'),'" "',...
                        pp.outfn{i}(1:end-4),'.raw"'];
        end
    end
end

