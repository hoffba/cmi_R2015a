% ImageClass function
% Set DCE perfusion fitting options
function setFitOpts(self,varargin)
% setFitOpts() provides a GUI to input values
% setFitOpts('Name',value,....)
% setFitOpts(fitOpts) where fitOpts is a structure of fitting options

fopts = self.model.fitOpts;

if (nargin>1)
    nin = length(varargin);
    if (round(nin/2) == nin/2)
        for i = 1:nin/2
            ind = 2*i - 1;
            if isfield(fopts,varargin{ind})
                fopts.(varargin{ind}) = varargin{ind+1};
            elseif strcmp(varargin{ind},'voxfit')
                self.voxfit = logical(varargin{ind+1});
            end
        end
        self.model.setFitOptions(fopts);
    end
else
    maxiter = self.model.fitOpts.MaxIter;
    maxfunevals = self.model.fitOpts.MaxFunEvals;
    tolx = self.model.fitOpts.TolX;
    tolfun = self.model.fitOpts.TolFun;
    dis = self.model.fitOpts.Display;

    prompt = {'MaxIter','MaxFunEvals','TolX','TolFun',...
        'Display (off,iter,final)','VoxFit (0 or 1)'};
    defAns = {num2str(maxiter),num2str(maxfunevals),num2str(tolx),...
        num2str(tolfun),dis,num2str(self.voxfit)};
    answer = inputdlg(prompt,'Perfusion Fit Options',1,defAns);
    if ~isempty(answer)
        % Check for valid MaxIter
        tval = str2double(answer{1});
        if tval>0
            maxiter = round(tval);
        end
        % Check for valid MaxFunEvals
        tval = str2double(answer{2});
        if tval>0
            maxfunevals = round(tval);
        end
        % Check for valid TolX
        tval = str2double(answer{3});
        if tval>0
            tolx = tval;
        end
        % Check for valid TolFun
        tval = str2double(answer{4});
        if tval>0
            tolfun = tval;
        end
        % Check for valid Display
        if any(strcmp(answer{5},{'off','iter','final'}))
            dis = answer{5};
        end
        % Create new options structure
        self.model.setFitOptions(optimset('MaxIter',maxiter,...
                                'MaxFunEvals',maxfunevals,...
                                'Display',dis,... % 'off' , 'on' , 'final'
                                'TolX',tolx,...
                                'TolFun',tolfun));
        % Set VoxFit
        self.voxfit = logical(str2double(answer{6}));
    end
end