function [prmvals,cmap,labels,ha] = prm_rgb(vals,options,dsopt,ha)
% Calculates original RGB PRM with 95% CI along the unit line 
% vals = [xvals,yvals]; (two columns)
if nargin && (size(vals,2)==2)
    if (nargin > 1) && (length(options)==5)
        thresh = abs(options(1));
        txmin = options(2);
        txmax = options(3);
        tymin = options(4);
        tymax = options(5);
    else
        thresh = 55; % Default value
        txmin = 0;
        txmax = 100000;
        tymin = 0;
        tymax = 100000;
    end
    if (nargin<3) || isempty(dsopt)
        dsopt = true; % default to showing the scatterplot
    end

    prmvals = zeros(size(vals,1),1);
    tdiff = diff(vals,1,2);
    prmvals((tdiff>thresh) & (vals(:,1)>txmin) & (vals(:,2)>tymin) ...
        & (vals(:,1)<txmax) & (vals(:,2)<tymax)) = 3;% Red
    prmvals((tdiff<-thresh) & (vals(:,1)>txmin) & (vals(:,2)>tymin) ...
        & (vals(:,1)<txmax) & (vals(:,2)<tymax)) = 1;% Blue
    prmvals((tdiff <= thresh) & (tdiff >= -thresh) ...
        & (vals(:,1)>txmin) & (vals(:,2)>tymin) & (vals(:,1)<txmax) & (vals(:,2)<tymax)) = 2;% Green
    cmap = [1 0 0;0 1 0;0 0 1]; % [R;G;B]
    labels = {'Red','Green','Blue'};
    
    % Now show scatterplot if desired and save stats
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
                'Color',cmap(i,:),'MarkerSize',4,'MarkerFaceColor',cmap(i,:));
        end
        % Plot region lines:
        tx = [0 50000];
        plot(ha,tx,tx,'-r',...
                tx,tx+thresh,'-k',...
                tx+thresh,tx,'-k');
        hold(ha,'off');
        %set(ha,'Color',[0.5 0.5 0.5])
        title(ha,char(str))
        set(get(ha,'XLabel'),'String','Pre-')
        set(get(ha,'YLabel'),'String','Post-')
        %axis(ha,'equal')
%         xlim(ha,[txmin txmax]);
%         ylim(ha,[tymin tymax]);
        xlim(ha,[1000 2500]);
        ylim(ha,[1000 2500]);
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





