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
addParameter(p,'outfn','',@(x)ischar(x)||iscellstr(x));
addParameter(p,'tp','',@ischar);
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

if ~isempty(pp.tp) && (~isempty(pp.in) || pp.jac || pp.jacmat || pp.def)
    outchk = length(pp.in)==length(pp.outfn);
    np = length(pp.in);
    str = cell(np,1);
    for i = 1:np
        str{i} = [fullfile(self.elxdir,'transformix'),...
               ' -out "',pp.odir,'"',...
               ' -tp "',pp.tp,'"'];
        if ~isempty(pp.in)
            str{i} = [str{i},' -in "',pp.in{i},'"'];
        end
        if pp.jac
            str{i} = [str{i},' -jac all'];
            pp.jac = false; % only need to perform once
        end
        if pp.def
            str{i} = [str{i},' -def all'];
            pp.def = false; % only need to perform once
        end
        if pp.jacmat
            str{i} = [str{i},' -jacmat all'];
            pp.jacmat = false; % only need to perform once
        end
        if ~isempty(pp.outfn{i}) && outchk
            if ispc
                cpexec = 'copy';
            else
                cpexec = 'cp';
            end
            str{i} = [str{i},self.sepstr,cpexec,' "',fullfile(pp.odir,'result.mhd'),'" "',...
                pp.outfn{i},'"',self.sepstr,cpexec,' "',fullfile(pp.odir,'result.raw'),...
                '" "',pp.outfn{i}(1:end-3),'raw";'];
        end
    end
end

