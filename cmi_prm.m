function [prmvals,cmap,labels,limits] = cmi_prm(type,xvals,yvals,options)
% Performs PRM analysis on input values
% opt defines the type of PRM analysis to perform
%   0: Normal RGB analysis
%   1: 4-color COPD Lung
%   2: 7-color CT Lung
%   3: 9-color Benjamin's Multi-Modality
% xvals & yvals are 1D vectors of image values
% options catches any optional input parameters (e.g. 95%CI)

prmvals = []; cmap = []; labels = {};
if (nargin > 2)
    if (any(type == 1:6) && ~isempty(xvals) && (size(xvals,2) == 1) ...
            && (length(xvals)==length(yvals)))
        % Switch options:
        %   1: RGB (w/ CI)
        %   2: 4-color COPD
        %   3: 7-color COPD
        %   4: Cardiac
        %   5: Fibrosis
        %   6: BOS-Rejection
        switch type
            case 1 % RGB
                if (nargin > 3)
                    thresh = options(1);
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
                prmvals = 0*xvals;
                tdiff = yvals - xvals;
                prmvals((tdiff>thresh) & (xvals>txmin) & (yvals>tymin) ...
                    & (xvals<txmax) & (yvals<tymax)) = 3;% Red
                prmvals((tdiff<-thresh) & (xvals>txmin) & (yvals>tymin) ...
                    & (xvals<txmax) & (yvals<tymax)) = 1;% Blue
                prmvals((yvals-xvals <= thresh) & (tdiff>=-thresh) ...
                    & (xvals>txmin) & (yvals>tymin) & (xvals<txmax) & (yvals<tymax)) = 2;% Green
                cmap = [1 0 0;0 1 0;0 0 1]; % [R;G;B]
                labels = {'Red','Green','Blue'};
                limits = [];
            case 2 % 4-color COPD Lung
                % options: [x-thresh; y-thresh]
                if (nargin > 3) && (length(options)==3)
                    xthresh = options(1);
                    ythresh = options(2);
                    maxthresh = options(3);
                else
                    xthresh = -856;
                    ythresh = -950;
                    maxthresh = -500;
                end
                inds = (xvals > maxthresh) | (yvals > maxthresh);
                prmvals = 0*xvals;
                prmvals((xvals >= xthresh) & (yvals >= ythresh) & ~inds) = 4;% Green
                prmvals((xvals < xthresh)  & (yvals >= ythresh) & ~inds) = 3;% Yellow
                prmvals((xvals < xthresh)  & (yvals < ythresh) & ~inds)  = 2;% Red
                prmvals((xvals >= xthresh) & (yvals < ythresh) & ~inds)  = 1;% White
                cmap = [0 1 0;1 1 0;1 0 0;1 1 1];
                labels = {'PRM-Norm','PRM-fSAD','PRM-Emph','White'};
                limits = [-1024 , maxthresh ; -1024 , maxthresh];
            case 3 % 7-color CT Lung
                % options: [CI; x-thresh; y-thresh]
                if (nargin == 4) && (length(options)==4)
                    CI = options(1);
                    xthresh = options(2);
                    ythresh = options(3);
                    maxthresh = options(4);
                else
                    CI = 94;
                    xthresh = -856;
                    ythresh = -950;
                    maxthresh = -500;
                end
                inds = (xvals >= maxthresh) | (yvals >= maxthresh);
                prmvals = 0*xvals;
                prmvals((yvals-xvals>=-CI) & (xvals>=xthresh) & (yvals>=ythresh) & ~inds) = 7;% Red
                prmvals((yvals-xvals<=-CI) & (xvals>=xthresh) & (yvals<ythresh) & ~inds) =  6;% Green
                prmvals((yvals-xvals<-CI)  & (xvals>=xthresh) & (yvals>=ythresh) & ~inds) = 5;% Blue
                prmvals((yvals-xvals>-CI)  & (xvals>=xthresh) & (yvals<ythresh) & ~inds) =  4;% Yellow
                prmvals((yvals-xvals<=-CI) & (xvals<xthresh)  & (yvals<ythresh) & ~inds) =  3;% Magenta
                prmvals((yvals-xvals>-CI)  & (xvals<xthresh)  & (yvals<ythresh) & ~inds) =  2;% Cyan
                prmvals((yvals-xvals>-CI)  & (xvals<xthresh)  & (yvals>=ythresh) & ~inds) = 1;% White
                cmap = [1 0 0;0 1 0;0 0 1;1 1 0;1 0 1;0 1 1;1 1 1];
                labels = {'Red','Green','Blue','Yellow','Magenta','Cyan','White'};
                limits = [-1024 , maxthresh ; -1024 , maxthresh];
            case 4 % Cardiac - shifted line RGB
                % options: [CI; max-thresh; slope; intercept]
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
                inds = (xvals<=thmax) & (yvals<=thmax) & (xvals>=thmin) & (yvals>=thmin);
                prmvals = 0*xvals;
                tdiff = yvals - slp*xvals - intcpt;
                prmvals((tdiff>CI) & inds) = 3;% Red
                prmvals((tdiff<-CI) & inds) = 1;% Blue
                prmvals((tdiff <= CI) & (tdiff>=-CI) & inds) = 2;% Green
                cmap = [1 0 0;0 1 0;0 0 1]; % [R;G;B]
                labels = {'Red','Green','Blue'};
                limits = [];
            case 5 % Fibrosis
                if (upythresh>lowythresh) && (upxthresh>lowxthresh)
                    nsub = 50000;
                    skip = (length(xvals)-1)/nsub;
                    xmin = -1024;
                    xmax = 200;
                    ymin = -1024;
                    ymax = 200;
                    % Make seperate plot for dscatter
                    figure;
                    h(1) = rectangle('Position',[xmin ymin (lowxthresh-xmin) (lowythresh-ymin)],'FaceColor','r');
                    h(2) = rectangle('Position',[xmin lowythresh (lowxthresh-xmin) (ymax-lowythresh)],'FaceColor','y');
                    h(3) = rectangle('Position',[lowxthresh lowythresh (xmax-lowxthresh) (ymax-lowythresh)],'FaceColor','g');
                    set(h,'LineWidth',2)
                    hold on
                    dscatter(xvals(round(1:skip:end)),yvals(round(1:skip:end)));
                    hold off
                    xlim([xmin xmax])
                    ylim([ymin ymax])
                    xlabel('Expiration (HU)')
                    ylabel('Inspiration (HU)')
                    %line([lowxthresh lowxthresh],[ymin ymax],'Color','k') % Lower X-threshold
                    h(1) = line([upxthresh upxthresh],[ymin ymax]); % Upper X-threshold
                    %line([xmin xmax],[lowythresh lowythresh],'Color','k') % Lower Y-threshold
                    h(2) = line([xmin xmax],[upythresh upythresh]); % Upper Y-threshold
                    h(3) = line([max(xmin,xmin-diag) min(xmax,xmax-diag)],...
                                [max(ymin,ymin+diag) min(ymax,ymax+diag)]); % Diagonal CI
                    set(h,'Color','k','LineWidth',2)
                end
            case 6 % BOS-Rejection
                % options: CI, exp/x thresh GT, insp/y thresh Emph, 
                % insp/y pfibr thresh -750,
                % parenchymal thresh/pre-fibrotic max cutoff -500,
                % fibrotic cutoff for ins/exp/x/y assuming vessels masked out 
                % Defaults: 94,-856,-950,-750,-500,0
                
                
                if (nargin > 3) && (length(options)==6)
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
                fprintf('cmi_prm: prmrejii prm values = %d\n',numel(xvals));
                
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
                
                prmvals = 0*xvals;
                % -- Overwrites certain values at each step to result in the final diagram above
                prmvals((xvals<=allmax)  & (yvals<=allmax))    = 11; % Dark Purple = visible fibrosis
                prmvals((xvals<=pthresh) & (yvals<=pthresh))   = 9;  % Dark Pink = Red COPD 1 but above F thresh
                prmvals((prmvals==9)     & (xvals<xlow))       = 10; % Light Pink = White COPD 3 but above F thresh
                prmvals((prmvals==9)     & (yvals<yhigh))      = 7;  % Orange
                prmvals((prmvals==10)    & (yvals<yhigh))      = 1;  % Yellow
                prmvals((prmvals==9)     & (yvals<(xvals-CI))) = 8;  % Light Purple
                prmvals((prmvals==7)     & (yvals<ylow))       = 6;  % White
                prmvals((prmvals==1)     & (yvals<ylow))       = 2;  % Red
                prmvals((prmvals==2)     & (yvals<(xvals-CI))) = 3;  % Red2
                prmvals((prmvals==7)     & (yvals<(xvals-CI))) = 5;  % Green
                prmvals(((prmvals==1) & (yvals<(xvals-CI))) ...
                      | ((prmvals==6) & (yvals>=(xvals-CI)))) = 4; % n/a (small area where xlow and ylow meet)
                
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
                limits = [-1024 , allmax ; -1024 , allmax];
        end
    end
end

