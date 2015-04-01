classdef HistogramClass < handle
    properties (SetAccess=private, GetAccess=public)
        active = true;      % Visible status of plot figure
        plotType = 0;       % Determination for type of histogram plot:
                            %   0 - Normal: Volume/Slice
                            %   1 - Multi-plot (e.g. timeseries)
                            %   2 - PRM
        nbins = 100;        % Number of histogram bins
        
        hImgFig = nan;      % Handle to image figure
        hfPlot = nan;       % Handle to plot figure
        haPlot = nan;       % Handle to plot axes
        plotFigPos          % position of plot figure
        
        hLines              % Handles of individual lineseries objects
        labels = {};        % Display labels for each line (cell array of strings)
        binVals             % Bin counts (matrix - column vectors)
        binLocs             % Locations of centers of bins (vector)
        
        % Statistics Tracking
        voxvol = 1;         % Image voxel volume
        means               % Vector of dataset Means
        stdevs              % Vector of dataset StdDevs
        meds                % Vector of dataset Medians
        
        % Display options
        lineStyle = '-';    % LineStyle for all line plots
        lineMark = 'none';  % Marker for line points
    end
    methods (Access = private)
        % Initialize plot figure
        function initPlotFig(self)
            if ~ishandle(self.hfPlot)
                self.hfPlot = figure('Position',self.plotFigPos);
                self.haPlot = axes('Parent',self.hfPlot,'YLim',[0 1]);
                hold all;
                title({[];[]});
                addlistener(self.hfPlot,'LocationChanged',@self.changeDispPos);
                self.dispUDvisible;
                if ~isempty(self.hLines)
                    self.hLines = plot(self.binLocs,self.binVals,...
                        'Parent',self.haPlot,'DisplayName',self.labels);
                end
            end
        end
        % Set visible state of plot figure based on active status
        function dispUDvisible(self)
            self.initPlotFig;
            if ishandle(self.hfPlot)
                if (self.active && ~isempty(self.hLines))
                    str = 'On';
                else
                    str = 'Off';
                end
                set(self.hfPlot,'Visible',str,'Position',self.plotFigPos);
            end
        end
        % Update printed statistics in title
        function dispUDstats(self,ind)
            % Optional input ind for Volume (1) or Slice (2)
            if (nargin == 1)
                ind = 0; % Update both
            end
            switch self.plotType
                case 0 % Normal: Volume/Slice
                    if (ind ~= 0)
                        str = get(get(self.haPlot,'Title'),'String');
                    else
                        str = {[];[]};
                    end
                    vols = self.voxvol * sum(self.binVals,1);
                    if any(ind == [0 1])
                        n = find(strcmp('Volume',self.labels),1);
                        if ~isempty(n)
                            str{1} = ['Volume: ' num2str(vols(n)) ...
                                   ' ; Mean= ' num2str(self.means(n)) ...
                                   ' ; StDev= ' num2str(self.stdevs(n)) ...
                                   ' ; Med= ' num2str(self.meds(n))];
                        end
                    end
                    if any(ind == [0 2])
                        n = find(strcmp('Slice',self.labels),1);
                        if ~isempty(n)
                            str{2} = ['Slice: ' num2str(vols(n)) ...
                                   ' ; Mean= ' num2str(self.means(n)) ...
                                   ' ; StDev= ' num2str(self.stdevs(n)) ...
                                   ' ; Med= ' num2str(self.meds(n))];
                        end
                    end
                    title(self.haPlot,str);
                case 1 % Multi-plot
                    
                case 2 % PRM
            end
        end
        
        % Listener function for updating the plot figure position
        function changeDispPos(self,~,~)
            self.plotFigPos = get(self.hfPlot,'Position');
        end
    end
    methods
        % Constructor with option to initialized with data
        function obj = HistogramClass(hImgFig,histMin,histMax,nbins,vals)
            % Determine plot figure position
            if nargin
                obj.hImgFig = hImgFig;
                str = get(hImgFig,'Units');
                if ~strcmp(str,'pixels')
                    set(hImgFig,'Units','pixels');
                end
                tpos = get(hImgFig,'Position');
                set(hImgFig,'Units',str);
                if (tpos(2) < 20) % Don't place plot under image if too low
                    obj.plotFigPos = get(0,'defaultfigureposition');
                else
                    % Default histogram plot figure height = 200 (pixels)
                    obj.plotFigPos = [tpos(1) tpos(2)-280 tpos(3) 200];
                end
            else
                obj.plotFigPos = get(0,'defaultfigureposition');
            end
            obj.initPlotFig;
            % Initialize histogram with optional inputs
            if (nargin == 3)
                obj.initialize(histMin,histMax);
            elseif (nargin == 4)
                obj.initialize(histMin,histMax,nbins);
            end
            if (nargin == 5) % Histogram the input values
                obj.binVals = hist(vals,obj.binLocs);
                self.hLines = plot(obj.binLocs,obj.binVals,...
                    'Parent',obj.haPlot,'DisplayName',obj.labels);
            end
        end
        % Destructor closes all related figures
        function delete(self)
            %disp('Deleting HistogramObject ... ')
            if ishandle(self.hfPlot)
                close(self.hfPlot)
            end
        end
        
    end
end