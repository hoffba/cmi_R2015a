% ElxClass function
% Generate system call to elastix
% Inputs: odir  = output directory
%         Name/Value pairs (need at least one of the following:)
%           'in'    / Images to transform (.mhd)
%           'jac'   / Flag to generate Jacobian map (TF)
%           'def'   / Flag to generate deformation maps (TF)
%           'jacmat'/ Flag to generate Jacobian matrices (TF)
%         Name/Value pairs (optional)
%           'outfn' / File names for transformed 'in' images
%           'tp'    / TransformParameterFile (.txt)
%               * if absent, uses ElxClass.Tx0 property
function str = tfxCmd(self,odir,varargin)

str = '';
% Parse variable inputs:
p = inputParser;
addRequired(p,'odir',@isfolder);
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
    % If no file was input, save current Tx0 in ElxClass object
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
nf = max(length(pp.in),(pp.jac || pp.jacmat || pp.def));
if length(pp.tp)~=nf
    pp.tp = repmat(pp.tp,nf,1);
end

if ~isempty(pp.tp) && (~isempty(pp.in) || pp.jac || pp.jacmat || pp.def)
    outchk = length(pp.outfn)==nf;
    str = cell(1+outchk,nf);
    for i = 1:nf
        str{1,i} = [fullfile(self.elxdir,'transformix'),...
               ' -out "',pp.odir,'"',...
               ' -tp "',pp.tp{i},'"'];
        if ~isempty(pp.in{i})
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
            outfn = pp.outfn{i};
            
            % Elastix can't save as zipped
            gzflag = endsWith(outfn,'.gz');
            if gzflag
                outfn = outfn(1:end-3);
            end

            [odir,outfn,ext] = fileparts(outfn);
            switch self.sys
                case 2 % PC
                    if strcmp(ext,'.mhd')
                        ext = {'.mhd','.raw'};
                    else 
                        ext = {ext};
                    end
                    for iext = 1:numel(ext)
                        str{2,i} = ['rename "',fullfile(odir,['result',ext{iext}]),'" "',outfn,ext{iext},'"'];
                    end
                case 4 % Great Lakes
                    str{2,i} = ['rename result ',outfn,' ',fullfile(odir,'result.???')];
                otherwise
                    str{2,i} = ['mv "',fullfile(odir,'result.???'),'" "',outfn,'.???"'];
            end
        end
    end
end
