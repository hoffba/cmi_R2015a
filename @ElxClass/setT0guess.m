% ElxClass function
function stat = setT0guess(self,str)

stat = false;
if (nargin==2) && isstruct(str) && all(isfield(str,{'i','fname','fpath'}))
    err = find(any(cellfun(@(x)~exist(x,'file'),cellfun(@(x,y)fullfile(x,y),...
        {str(:).fpath},{str(:).fname},'UniformOutput',false))),1);
    if err
        error(['File does not exist: ',str(err).fname]);
    end
    self.Tx0guess = str;
    stat = true;
end