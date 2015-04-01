% ModelClass function
% Set fitting options
function setFitOptions(self,options)
ask = true;
prompt = {'MaxIter','MaxFunEvals','TolX','TolFun','Display'};
if (nargin==2) && isstruct(options) && all(isfield(options,prompt))
    self.fitOpts = options;
    ask = false;
end
if ask
    maxiter = self.fitOpts.MaxIter;
    maxfunevals = self.fitOpts.MaxFunEvals;
    tolx = self.fitOpts.TolX;
    tolfun = self.fitOpts.TolFun;
    dis = self.fitOpts.Display;
    defAns = {num2str(maxiter),num2str(maxfunevals),num2str(tolx),...
        num2str(tolfun),dis};
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
        self.fitOpts = optimset('MaxIter',maxiter,...
                                'MaxFunEvals',maxfunevals,...
                                'Display',dis,... % 'off' , 'on' , 'final'
                                'TolX',tolx,...
                                'TolFun',tolfun);
    end
end