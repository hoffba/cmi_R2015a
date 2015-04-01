% ModelClass function
% Set initial guess
function setPar0(self,val)
tdef = self.defs.(self.getModType)(self.cmodel(self.mtype)).par0;
np = length(tdef);
if (nargin==2) && isnumeric(val) && (length(val(:))==np)
    self.defs.(self.getModType)(self.cmodel(self.mtype)).par0 = val(:)';
else
    def = cellfun(@num2str,mat2cell(tdef',ones(1,np))','UniformOutput',0);
    answers = str2double(inputdlg(self.defs.(self.getModType)(self.cmodel(self.mtype)).labels,...
        ['Initial Guess: (' self.getModName ')'],1,def));
    if ~isempty(answers) && ~any(isnan(answers))
        self.defs.(self.getModType)(self.cmodel(self.mtype)).par0 = answers(:)';
    end
end