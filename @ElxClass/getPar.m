% ElxClass function
% Set Elastix parameters in existing schedule
function varargout = getPar(self,i,varargin)
% Get current parameter value fldn from schedule step i

na = length(varargin);
if nargin<3
    error('Not enough input arguments.')
elseif ~ismember(i,1:length(self.Schedule))
    error('Invalid Schedule index.')
elseif ~iscellstr(varargin)
    error('Invalid parameter inputs.')
elseif nargout<na
    warning(['Not enough outputs for requested parameters ... only returning the first ',...
        num2str(nargout),'.']);
    varargin(nargout+1:end) = [];
elseif nargout>na
    warning(['Too many outputs ... only return the first ',num2str(na),'.']);
end
varargout = cell(1,nargout);
for j = find(isfield(self.Schedule{i},varargin))
    varargout{j} = self.Schedule{i}.(varargin{j});
end
