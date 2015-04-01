% Class containing methods for fitting data to a model
classdef FitClass < handle
    %FitClass Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=protected, GetAccess=public)
        typeStr = 'Exp';    % Type of model
        xdata               % [n x 1] Independent variable (time, b-value, etc.)
        ydata               % [n x 1] Dependent variable (SI, etc.)
        yfit                % Modeled data
        
        mind = 1;           % Index of selected model
        defs = struct('name','Exponential',...
                      'par0',[1,-0.5],...
                      'func',@self.fitExp,...
                      'labels',{'A','k'},...
                      'bounds',[0,inf;-5,5]);
        
        OptimChk = false;   % Check if Optimization toolbox is available for fits
        fitOpts = struct('MaxIter',    100,...
                         'MaxFunEvals', 1000,...
                         'Display',     'iter',...
                         'TolX',        1e-4,...
                         'TolFun',      1e-4,...
                         'VoxChk',      false);
        
        % Structure of handles to plot objects
        h = struct('fig',[],'ax',[],'dline',[],'fline',[],'title',[]);
        
    end
    methods (Access=private)
        % Minimization function for curve fitting
        function chi = fitFunc(self,par)
            if self.mind
                self.calcYFit(par);
                if strcmp(self.fitOpts.Display,'iter')
                    self.updatePlot('fit',self.yfit);
                end
                pause(0.005)
                % For use with LSQNONLIN:
                chi = (self.yfit - self.ydata).^2;
                % For use with FMINSEARCH:
                %chi = sum(chi);
            end
        end
        % Update plot values
        function updatePlot(self,str)
            go = true;
            switch str
                case 'data'
                    th = self.h.dline;
                    y = self.ydata;
                case 'fit'
                    th = self.h.fline;
                otherwise
            end
            if strcmp(str,'data')
            elseif strcmp(str,'fit') && (nargin==3)
            else
                go = false;
            end
            if go
                set(th,'YData',y);
            end
        end
    end
    methods
        % CONSTRUCTOR
        function self = FitClass
            % Check if Optimization Toolbox is available for fits
            [TF,~] = license('checkout','Optimization_Toolbox');
            if TF
                self.OptimChk = true;
            else
                warning('Curve-Fitting not available - unable to checkout Optimization Toolbox license!')
            end
        end
        % DESTRUCTOR
        function delete(self)
            if ishandle(self.h.fig)
                delete(self.h.fig)
            end
        end
        % Set xdata
        function stat = setXdata(self,x)
            stat = false;
            if nargin==1 % use GUI to input x-values
                x = [];
                answer = inputdlg('Input X-Data:','X-Data',1);
                if ~isempty(answer)
                    x = eval(answer);
                end
            end
            nx = length(x);
            if (nx>1) && isnumeric(x) && ~any(isnan(x))
                self.xdata = x;
                stat = true;
                % Make sure ydata is same length
                if length(self.ydata) ~= nx
                    self.ydata = zeros(1,nx);
                end
            end
        end
        % Set ydata
        function stat = setYData(self,y)
            stat = false;
            if (nargin==2) && ~isempty(y)
                self.ydata = y(:)';
                stat = true;
                if length(self.xdata) ~= length(self.ydata)
                    self.xdata = 0:(length(self.ydata)-1);
                end
            end
        end
        % Set current model
        function ok = setModel(self,val)
            if (nargin==1) && (self.mind>0)
                % Put list of models together for user selection
                opts = {self.defs(:).name};
                [val,ok] = listdlg('ListString',opts,'InitialValue',self.mind);
            end
            if isstruct(val) && all(isfield(val,{'name','par0','func','labels','bounds'})) ...
                    && isa(val.func,'function_handle')
                % Add to defined models
                self.defs(end+1) = val;
                self.mind = length(self.defs);
            elseif isnumeric(val) && any(val==(1:length(self.defs)))
                self.mind = val;
            else
                warning(['Invalid input for model type. Must be integer from 1 to ',...
                         num2str(length(self.defs))])
            end
        end
        % Set curve-fitting options
        function setFitOptions(self,varargin)
        % setFitOptions() uses GUI for inputs
        % setFitOptions('Name',value,...) sets options ('Name') to values
        % setFitOptions(optionsStruct) explicitly sets options
            prompt = {'MaxIter','MaxFunEvals','TolX','TolFun','Display','VoxChk'};
            if nargin==1
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
                        self.fitOpts.MaxIter = round(tval);
                    end
                    % Check for valid MaxFunEvals
                    tval = str2double(answer{2});
                    if tval>0
                        self.fitOpts.MaxFunEvals = round(tval);
                    end
                    % Check for valid TolX
                    tval = str2double(answer{3});
                    if tval>0
                        self.fitOpts.TolX = tval;
                    end
                    % Check for valid TolFun
                    tval = str2double(answer{4});
                    if tval>0
                        self.fitOpts.TolFun = tval;
                    end
                    % Check for valid Display
                    if any(strcmp(answer{5},{'off','iter','final'}))
                        self.fitOpts.Display = answer{5};
                    end
                    % Check for valid VoxChk value
                    tval = str2double(answer{6});
                    if ismember(tval,[0,1])
                        self.fitOpts.VoxChk = logical(tval);
                    end
                end
            elseif isstruct(varargin{1})
                for i = 1:length(prompt)
                    if isfield(varargin{1},prompt{i})
                        self.fitOpts.(prompt{i}) = varargin{1}.(prompt{i});
                    end
                end
            elseif nargin>2
                fn = varargin(1:2:end);
                vals = varargin(2:2:end);
                for i = 1:length(vals)
                    if ischar(fn{i}) && isfield(self.fitOpts,fn{i})
                        self.fitOpts.(fn{i}) = vals{i};
                    end
                end
            end
        end
        % Get current model's name
        function str = getModName(self)
            str = self.defs(self.mind).name;
        end
        % Generate current model's YData
        function y = calcYFit(self,par)
            if nargin==1
                par = self.defs(self.mind).par0;
            end
            y = feval(self.defs(self.mind).func,par);
            self.yfit = y;
        end
        % Calculate goodness of fit for current model
        function gof = calcGoF(self)
            if nargin==2
                gof = sqrt(sum((self.yfit-self.ydata).^2)/length(y)); % RMSE
                gof = gof/(max(self.ydata) - min(self.ydata)); % NRMSE
            end
        end
    end
end