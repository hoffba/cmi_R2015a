function [prmvals,cmap,labels,ha] = prm_cardiac(vals,options,dsopt,ha)
% Calculates shifted-line PRM for Cardiac
% vals = [xvals,yvals]; (two columns)
% options: [CI; max-thresh; slope; intercept]
if (nargin>0) && (size(vals,2)==2)
    if (nargin == 4) && (length(options)==5)
        CI = options(1);
        thmax = options(2);
        thmin = options(3);
        slp = options(4);
        intcpt = options(5);
    else
        CI = 2000;
        thmax = 20000;
        thmin = 500;
        slp = 0.5;
        intcpt = 2200;
    end
    if (nargin<3) || isempty(dsopt)
        dsopt = true; % default to showing the scatterplot
    end
    inds = (vals<=thmax) & (vals>=thmin);
    inds = inds(:,1) & inds(:,2);
    
    prmvals = zeros(size(vals,1),1);
    tdiff = vals(:,2) - slp*vals(:,1) - intcpt;
    prmvals((tdiff>CI) & inds) = 3;% Red
    prmvals((tdiff<-CI) & inds) = 1;% Blue
    prmvals((tdiff <= CI) & (tdiff>=-CI) & inds) = 2;% Green
    cmap = [1 0 0;0 1 0;0 0 1]; % [R;G;B]
    labels = {'Red','Green','Blue'};
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
        hold(ha,'off');
        %set(ha,'Color',[0.5 0.5 0.5])
        title(ha,char(str))
        set(get(ha,'XLabel'),'String','Diastole')
        set(get(ha,'YLabel'),'String','Systole')
        %xlim(ha,[600 1800]);
        %ylim(ha,[600 1800]);
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