function [prmvals,cmap,labels,ha] = prm_Dscat(vals,options,dsopt,ha)
% Generates a density scatterplot
% Inputs:
%       vals:       [xvals,yvals]; (two columns)
%       options:    [axis-min; axis-max]
%       dsopt:      check to display the scatterplot
%       ha:         handle of existing scatterplot
% Outputs:
%       prmvals:    values of PRM voxels
%       cmap:       colormaps corresponding to PRM regions
%       labels:     labels corresponding PRM regions
%       ha:         handle of scatterplot

if (nargin>0) && (size(vals,2)==2)
    if (nargin > 3) && (length(options)==3)
        alim = [options(1),options(2)];
        dscale = options(3);
    else
        alim = [-1024,200];
        dscale = nan;
    end
    if (nargin<3) || isempty(dsopt)
        dsopt = true; % default to showing the scatterplot
    end
    
    np = length(vals);
    [prmvals,sfact] = calcDscatter(vals(:,1),vals(:,2)); % Calculate on full data set
    if dscale==0
        dscale = sfact;
    else
        prmvals = prmvals * sfact / dscale;
    end
    cmap = jet(64);
    labels = {};
    if dsopt
        if (nargin<4) || isempty(ha) || ~(ishandle(ha) && strcmp(get(ha,'Tag'),'PRMscatter'))
            figure;
            ha = axes('Tag','PRMscatter');
        else
            cla(ha);
        end
        np2 = 10000; % number of desired points to scatterplot
        if np > np2
            % sub-sampled PRM indexes (target 10k points)
            tp = unique(round(1:(np/np2):np));
        else
            tp = 1:np;
        end
        xvals = vals(tp,1);
        yvals = vals(tp,2);
        D = prmvals(tp,1);
        scatter(ha,xvals,yvals,20,D,'filled')
        if ~isnan(dscale)
            caxis(ha,[0,1])
        end
        hold(ha,'on'), plot(ha,alim,alim,'-k'), hold(ha,'off')
        set(ha,'Tag','PRMscatter') % scatter function must clear this property
        title(ha,'Density Scatterplot')
        set(get(ha,'XLabel'),'String','Expiration')
        set(get(ha,'YLabel'),'String','Inspiration')
        xlim(ha,alim);
        ylim(ha,alim);
    end
    prmvals = round(prmvals*64);
end