% ModelClass function
% Allows modification of DCE-specific parameters
function setDCEopts(self,varargin)
% INPUT: varargin = (optional) name/value pairs for valid options

ok = true;
if (nargin==1) % Use GUI to input values
    fnames = fieldnames(self.dceOpts);
    tdefs = cellfun(@num2str,struct2cell(self.dceOpts),'UniformOutput',0);
    % Display input dialog box
    vals = cellfun(@str2num,inputdlg(fnames,['Options - ',self.getModName],1,tdefs));
elseif isstruct(varargin{1})
    fnames = fieldnames(varargin{1});
    vals = cell2mat(struct2cell(varargin{1}));
elseif nargin>2 % Parse values for valid input options
    fnames = varargin(1:2:end);
    vals = varargin{2:2:end};
else
    ok = false;
end
if ok && (length(fnames)==length(vals))
    for i = 1:length(fnames)
        if isfield(self.dceOpts,fnames{i}) && isnumeric(vals(i)) && ~isnan(vals(i))
            self.dceOpts.(fnames{i}) = vals(i);
            if strcmp(fnames{i},'flipa')
                % Also set cosine and sine
                self.cf = cosd(self.dceOpts.flipa);
                self.sf = sind(self.dceOpts.flipa);
            end
        else
            warning(['Invalid input: ',char(fnames{i})]);
        end
    end
end