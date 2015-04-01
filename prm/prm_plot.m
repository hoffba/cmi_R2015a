function [ha,hs] = prm_plot(x,y,ind,cmap,labels,pcts,ha,varargin)
% Inputs: x      = x-coordinates
%         y      = y-coordinates
%         ind    = prm index
%         cmap   = index colors
%         labels = cell array of string labels
%         pcts   = matrix of region percents
%         ha     = handle of existing scatterplot axes
%         varargin = optional inputs: Namve/Value pairs

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
        p.addParamValue('size',      12,@(x) isnumeric(x) && (x>0));
        p.addParamValue('fill',      true,@(x) islogical(x));
        p.addParamValue('markertype','o',@(x) ischar(x));
        p.addParamValue('xlim',      'auto',@(x) isempty(x) || ...
                                     (isnumeric(x) && all(size(x)==[1,2])));
        p.addParamValue('ylim',      'auto',@(x) isempty(x) || ...
                                     (isnumeric(x) && all(size(x)==[1,2])));
        p.addParamValue('xlabel',    '',@ischar);
        p.addParamValue('ylabel',    '',@ischar);
        p.addParamValue('bgcolor',   0.6*[1,1,1],@(x) isnumeric(x) ...
                                                 && all(size(x)==[1,3]) ...
                                                 && all(x>=0) && all(x<=1));
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
        colormap(ha,cmap);
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
               'PlotBoxAspectRatioMode','manual');
        set(get(ha,'XLabel'),'String',opts.xlabel);
        set(get(ha,'YLabel'),'String',opts.ylabel);
        xlim(ha,opts.xlim);
        ylim(ha,opts.ylim);
        
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
