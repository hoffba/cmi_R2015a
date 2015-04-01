function [prmvals,cmap,labels,ha] = prm_copd4(vals,options,dsopt,ha)
% Generates a density scatterplot
% Inputs:
%       vals:       [xvals,yvals]; (two columns)
%       options:    [axis-min; axis-max]
%       dsopt:      check to display the scatterplot
%       ha:         handle of existing scatterplot
if (nargin>0) && (size(vals,2)==2)
    if (nargin > 3) && (length(options)==2)
        xthresh = options(1);
        ythresh = options(2);
        maxthresh = options(3);
    else
        xthresh = -856;
        ythresh = -950;
        maxthresh = -500;
    end
    if (nargin<3) || isempty(dsopt)
        dsopt = true; % default to showing the scatterplot
    end
    inds = (vals >= maxthresh); inds = (inds(:,1) | inds(:,2));
    prmvals = zeros(size(vals,1),1);
    prmvals((vals(:,1) >= xthresh) & (vals(:,2) >= ythresh) & ~inds) = 4;% Green
    prmvals((vals(:,1) < xthresh)  & (vals(:,2) >= ythresh) & ~inds) = 3;% Yellow
    prmvals((vals(:,1) < xthresh)  & (vals(:,2) < ythresh) & ~inds)  = 2;% Red
    prmvals((vals(:,1) >= xthresh) & (vals(:,2) < ythresh) & ~inds)  = 1;% White
    cmap = [0 1 0;1 1 0;1 0 0;1 1 1];
    labels = {'PRM-Norm','PRM-fSAD','PRM-Emph','White'};
    % Now show scatterplot if desired
    np = sum(prmvals(:)>0); % only use non-zero values for total volume
    nprm = length(labels);
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
        cvals = prmvals(tp);
        str = cell(1,nprm);
        prmres = zeros(1,nprm);
        hold(ha,'on');
        for i = 1:nprm
            ii = nprm-i+1;
            prmres(i) = sum(prmvals(:)==ii) / np * 100;
            str{i} = [labels{i} '= ' num2str(prmres(i)) ' ; '];
            plot(ha,xvals(cvals == ii),yvals(cvals == ii),'o',...
                'Color',cmap(i,:),'MarkerSize',2,'MarkerFaceColor',cmap(i,:));
        end
        % Plot region lines:
        plot(ha,[xthresh xthresh],[-1024 maxthresh],'-k',...
                [-1024 maxthresh],[ythresh ythresh],'-k');
        hold(ha,'off');
        set(ha,'Color',0.6*[1,1,1])
        title(ha,char(str))
        set(get(ha,'XLabel'),'String','Expiration')
        set(get(ha,'YLabel'),'String','Inspiration')
        xlim(ha,[-1024 maxthresh]);
        ylim(ha,[-1024 maxthresh]);
    else
        prmres = zeros(1,nprm);
        for i = 1:nprm
            prmres(i) = sum(prmvals(:)==i) / np * 100;
        end
        prmres = fliplr(prmres);
    end
    assignin('base','prmpercents',prmres);
    assignin('base','prmlabels',labels);
    assignin('base','prmvals',[xvals yvals]);
end