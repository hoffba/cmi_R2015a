function [ha,hs] = prm_plot(x,y,ind,cmap,labels,pcts,ha,varargin)
% Inputs: x      = x-coordinates
%         y      = y-coordinates
%         ind    = prm index
%         cmap   = index colors
%         labels = cell array of string labels
%         pcts   = matrix of region percents
%         ha     = handle of existing scatterplot axes
%         optional inputs (Name/Value pairs):
%                   'size'          marker size
%                   'fill'          marker fill
%                   'markertype'    marker type
%                   'xlim'          [1x2] x-axis limits
%                   'ylim'          [1x2] y-axis limits
%                   'xlabel'        x-axis label
%                   'ylabel'        y-axis label
%                   'bgcolor'       axes background color [R G B]
%                   'm'             [1xn] threshold line slopes
%                   'b'             [1xn] threshold line y-intercepts

if (nargin>=4)
    nx = length(x);
    ny = length(y);
    ni = length(ind);
    nc = size(cmap,1);
    if (nx>0) && (nx==ny) && (ny==ni)
        % If no labels input, create some
        if (nargin<5) || isempty(labels)
            labels = strcat({'Area'},strtrim(cellstr(num2str(1:nc)')))';
        end
        % If no stats input, generate them
        % (coord vectors may have been sub-sampled)
        if (nargin<6) || isempty(pcts)
            for i = 1:nc
                pcts(i) = sum(ind==i) / nx * 100;
            end
        end
        % If no axes handle was input, create a new figure
        if (nargin<7) || isempty(ha) || ~(ishandle(ha) && ...
                    strcmp(get(ha,'Tag'),'PRMscatter'))
            figure('Name','PRMscatter'), ha = axes;
        else % Clear current axes
            cla(ha);
        end
        p = inputParser;
        p.addParameter('size',      12,@(x) isnumeric(x) && (x>0));
        p.addParameter('fill',      true,@(x) islogical(x));
        p.addParameter('markertype','o',@(x) ischar(x));
        p.addParameter('xlim',      'auto',@(x) isempty(x) || ...
                                     (isnumeric(x) && all(size(x)==[1,2])));
        p.addParameter('ylim',      'auto',@(x) isempty(x) || ...
                                     (isnumeric(x) && all(size(x)==[1,2])));
        p.addParameter('xlabel',    '',@ischar);
        p.addParameter('ylabel',    '',@ischar);
        p.addParameter('bgcolor',   0.6*[1,1,1],@(x) isnumeric(x) ...
                                                 && all(size(x)==[1,3]) ...
                                                 && all(x>=0) && all(x<=1));
        p.addParameter('m',[],@isvector);
        p.addParameter('b',[],@isvector);
        if (nargin>8)
            p.parse(varargin{:});
        end
        opts = p.Results;
        if isempty(opts.xlim)
            opts.xlim = 'auto';
        end
        if isempty(opts.ylim)
            opts.ylim = 'auto';
        end
        
        % Remove cropped values (ind==0)
        x(~ind) = [];
        y(~ind) = [];
        ind(~ind) = [];
        
        % ~~~~~~ Individual point colors ~~~~~~~~~
        % Determine voxel colors:
%         C = zeros(length(x),3);
%         for i = 1:size(cmap,1)
%             ii = ind==i;
%             C(ii,:) = repmat(cmap(i,:),sum(ii),1);
%         end
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % ~~~~~~ Colormap color control ~~~~~~~~~
        colormap(ha.Parent,cmap);
        C = ind;
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Generate scatter plot
        if opts.fill
            hs = scatter(ha,x,y,opts.size,C,opts.markertype,'fill');
        else
            hs = scatter(ha,x,y,opts.size,C,opts.markertype);
        end
        set(hs,'Tag','PRMscatter');
            
        % Adjust axes properties
        set(ha,'Color',opts.bgcolor,'Tag','PRMscatter',...
               'CLim',[0 , size(cmap,1)]+0.5,...
               'PlotBoxAspectRatioMode','manual','FontSize',18);
        set(get(ha,'XLabel'),'String',opts.xlabel,'Interpreter','none',...
            'FontSize',22,'FontWeight','bold');
        set(get(ha,'YLabel'),'String',opts.ylabel,'Interpreter','none',...
            'FontSize',22,'FontWeight','bold');
        
        % Axes limits:
        if ischar(opts.xlim) || isinf(opts.xlim(1)) || isnan(opts.xlim(1))
            tlim(1) = min(x);
        else
            tlim(1) = opts.xlim(1);
        end
        if ischar(opts.xlim) || isinf(opts.xlim(2)) || isnan(opts.xlim(2))
            tlim(2) = max(x);
        else
            tlim(2) = opts.xlim(2);
        end
        xlim(ha,tlim);
        if ischar(opts.ylim) || isinf(opts.ylim(1)) || isnan(opts.ylim(1))
            tlim(1) = min(y);
        else
            tlim(1) = opts.ylim(1);
        end
        if ischar(opts.ylim) || isinf(opts.ylim(2)) || isnan(opts.ylim(2))
            tlim(2) = max(y);
        else
            tlim(2) = opts.ylim(2);
        end
        ylim(ha,tlim);
        
        % Plot the threshold lines:
        if ~isempty(opts.m) && (length(opts.m)==length(opts.b))
            nl = length(opts.m);
            plotlines(ha,[opts.m,1],[opts.b,0],[repmat({'-k'},1,nl),{'--k'}]);
        end
        
%         % Generate title w/ stats
%         nstr = size(char(labels(:)),2);
%         str = char(cellfun(@(str,val) ...
%                                 sprintf(['%-',num2str(nstr),'s = % 7.3f'],...
%                                         str,val),...
%                            labels(:),num2cell(pcts(:)),'UniformOutput',false));
%         title(ha,char(str),'FontName','FixedWidth');
    end
else
    ha = []; hs = [];
end
