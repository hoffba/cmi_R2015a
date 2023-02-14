function setOpts(self,varargin)

p = inputParser;
addParameter(p,'ext',           nan,    @ischar);
addParameter(p,'gl_flag',       nan,    @islogical);
addParameter(p,'username',      nan,    @ischar);
addParameter(p,'orient_flag',   nan,    @islogical);
addParameter(p,'swap_flag',     nan,    @islogical);
addParameter(p,'seg_method',    nan,    @islogical);
addParameter(p,'unreg_flag',    nan,    @islogical);
addParameter(p,'airway_flag',   nan,    @islogical);
addParameter(p,'scatnet_flag',  nan,    @islogical);
addParameter(p,'vessel_flag',   nan,    @islogical);
addParameter(p,'reg_flag',      nan,    @islogical);
addParameter(p,'quickreg',      nan,    @islogical);
addParameter(p,'prm_flag',      nan,    @islogical);
addParameter(p,'tprm_flag',     nan,    @islogical);
addParameter(p,'yacta',         nan,    @iscellstr);
addParameter(p,'prm',           nan,    @(x)iscell(x)||isstruct(x));
addParameter(p,'tprm',          nan,    @(x)iscell(x)||isstruct(x));
parse(p,varargin{:});

fldnames = fieldnames(p.Results);
for i = 1:numel(fldnames)
    tfld = fldnames{i};
    if ~isnan(p.Results.(tfld))
        switch tfld
            case 'yacta'
            case 'prm'
            case 'tprm'
            otherwise
                self.(tfld) = p.Results.(tfld);
        end
    end
end
