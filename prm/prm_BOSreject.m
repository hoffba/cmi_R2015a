function [prmvals,cmap,labels,ha] = prm_BOSreject(vals,options,dsopt,ha)
% Calculates BOS-Rejection PRM for detecting fibrosis
% vals = [xvals,yvals]; (two columns)
% options: CI, exp/x thresh GT, insp/y thresh Emph, 
% insp/y pfibr thresh -750,
% parenchymal thresh/pre-fibrotic max cutoff -500,
% fibrotic cutoff for ins/exp/x/y assuming vessels masked out 
% Defaults: 94,-856,-950,-750,-500,0

if (nargin>0) && (size(vals,2)==2)
    if (nargin==4) && (length(options)==6)
        CI = abs(options(1));
        xlow = options(2);
        ylow = options(3);
        yhigh = options(4); % pre-fibrotic above this val
        pthresh = options(5); % parenchymal maximum
        allmax = options(6); % fibrosis maximum
    else
        CI = 94;
        xlow = -856;
        ylow = -950;
        yhigh = -750;
        pthresh = -500;
        allmax = 200;
    end
    if (nargin<3) || isempty(dsopt)
        dsopt = true; % default to showing the scatterplot
    end
    fprintf('cmi_prm: prmrejii prm values = %d\n',length(vals));

    %|______________________________________________|___
    %|                                              |
    %|                                              |
    %|                                              |
    %|                                              |
    %|______________________                        |
    %|     |               |                        |
    %|  10 |    9         /|                        |
    %|     |            /  |           11           |
    %|_____|__________/__8_|                        |
    %|     |        /      |                        |
    %|  1  |  7   /        |                        |
    %|     |    /     5    |                        |
    %|     |  /            |                        |
    %|_____|/______________|                        |
    %| 2  /|       6       |                        |
    %|__/_3|_______________|________________________|__

    tdiff = diff(vals,1,2);
    prmvals = zeros(size(vals,1),1);
    % -- Overwrites certain values at each step to result in the final diagram above
    prmvals((vals(:,1)<=allmax)  & (vals(:,2)<=allmax))    = 11; % Dark Purple = visible fibrosis
    prmvals((vals(:,1)<=pthresh) & (vals(:,2)<=pthresh))   = 9;  % Dark Pink = Red COPD 1 but above F thresh
    prmvals((prmvals==9)     & (vals(:,1)<xlow)) = 10;  % Light Pink = White COPD 3 but above F thresh
    prmvals((prmvals==9)     & (vals(:,2)<yhigh)) = 7;  % Orange
    prmvals((prmvals==10)    & (vals(:,2)<yhigh)) = 1;  % Yellow
    prmvals((prmvals==9)     & (tdiff<-CI))        = 8;  % Light Purple
    prmvals((prmvals==7)     & (vals(:,2)<ylow))  = 6;  % White
    prmvals((prmvals==1)     & (vals(:,2)<ylow))  = 2;  % Red
    prmvals((prmvals==2)     & (tdiff<-CI))        = 3;  % Red2
    prmvals((prmvals==7)     & (tdiff<-CI))         = 5;  % Green
    prmvals(((prmvals==1) & (tdiff<-CI)) ...
          | ((prmvals==6) & (tdiff>=-CI))) = 4; % n/a (small area where xlow and ylow meet)

    lungvolume = sum(prmvals~=0);         
    fprintf('cmi_prm: Lung volume without fibrosis prmrejii (pixels) = %d\n',lungvolume);

    cmap = [0.4 , 0 , 0.7;...
            1 , 0.5 , 1;...
            1 , 0 , 1;...
            0.7 , 0 , 0.7;...
            1 , 0.5 , 0;...
            1 , 1 , 1;...
            0 , 1 , 0;...
            0 , 0 , 0;...
            1 , 0 , 0;...
            1 , 0 , 0;...
            1 , 1 , 0];
    labels = {'PRM Fibrosis 11',...
              'PRM EarlyFibrosis 10',...
              'PRM EarlyFibrosis 9',...
              'PRM EarlyFibrosis 8',...
              'PRM Normal 2 subF',...
              'PRM 6 ignored',...
              'PRM Normal 1 subF',...
              'PRM null 7',...
              'PRM Emph 5',...
              'PRM Emph 4',...
              'PRM fSAD 3 subF'}; 
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
        % ~~ Density scatter settings
        D = calcDscatter(xvals,yvals);
        nshades = 11; shlim = [1,0.4]; dsh = (shlim(1)-shlim(2))/(nshades-1);
        dcmap = zeros(nshades*nprm,3);
        tpad = 0.05;
        % ~~~ end ~~~
        hold(ha,'on');
        for i = 1:nprm
            ii = nprm-i+1;
            prmres(i) = sum(prmvals(:)==ii) / np * 100;
            str{i} = [labels{i} '= ' num2str(prmres(i)) ' ; '];
            ti = nshades*(i-1) + 1;
            % ~~~ Density colormap construction:
            dcmap(ti:(ti+nshades-1),:) = (shlim(1):-dsh:shlim(2))'*cmap(i,:);
            % ~~~ Plot straight colors:
            % plot(ha,xvals(cvals == ii),yvals(cvals == ii),'o',...
            %     'Color',cmap(i,:),'MarkerSize',2,'MarkerFaceColor',cmap(i,:));
        end
        % ~~~ Plot the density color-adjusted points:
        scatter(ha,xvals,yvals,20,(cvals-0.5+(1-2*tpad)*(D.^2-0.5))/(nprm),'filled')
        colormap(ha,flipud(dcmap));
        set(ha,'CLim',[0 1]);
        % ~~~~~~ end ~~~~~~~~
        % Plot region lines:
        plot(ha,[CI-1024,pthresh],[-1024,pthresh-CI],'-k',...
                [xlow,xlow],[-1024,pthresh],'-k',...
                [pthresh,pthresh],[-1024,pthresh],'-k',...
                [-1024,pthresh],[pthresh,pthresh],'-k',...
                [-1024,pthresh],[ylow,ylow],'-k',...
                [-1024,pthresh],[yhigh,yhigh],'-k');
        hold(ha,'off');
        set(ha,'Color',0.8*[1,1,1]) % gray BG to better see the dots
        title(ha,char(str))
        set(get(ha,'XLabel'),'String','Expiration')
        set(get(ha,'YLabel'),'String','Inspiration')
        xlim(ha,[-1024 allmax]);
        ylim(ha,[-1024 allmax]);
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
