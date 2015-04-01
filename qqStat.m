function stats = qqStat(varargin)

fnames = dir('*.mat');
if length(fnames)~=2
    fnames = uigetfile('*.mat','Select MAT files with PRM data:','',...
        'Multiselect','on');
else
    fnames = struct2cell(fnames);
    fnames = fnames(1,:);
end
if ~isempty(fnames) && iscell(fnames)
    if length(fnames)==2
        good = true;
        prmvals = open(fnames{1});
        if all(isfield(prmvals,{'xvals','yvals'}))
            x1 = prmvals.xvals;
            y1 = prmvals.yvals;
        else
            good = false;
        end
        if good
            prmvals = open(fnames{2});
            if all(isfield(prmvals,{'xvals','yvals'}))
                x2 = prmvals.xvals;
                y2 = prmvals.yvals;
            else
                good = false;
            end
            clear prmvals
            if good
                alims = [-1000 0];
                stats = zeros(4,4);
                figure;
                % First PRM, Ins/Exp
                subplot(2,2,1);
                [~,~,stats(1,:)] = getQQdata(x1,y1);
                xlim(alims),ylim(alims)
                title('X1-Y1','FontWeight','bold','FontSize',12)
                text(-900,-100,['R^2 = ' num2str(stats(1,1))])
                text(-900,-200,['p = ' num2str(stats(1,3))])
                % Second PRM, Ins/Exp
                subplot(2,2,2)
                [~,~,stats(2,:)] = getQQdata(x2,y2);
                xlim(alims),ylim(alims)
                title('X2-Y2','FontWeight','bold','FontSize',12)
                text(-900,-100,['R^2 = ' num2str(stats(2,1))])
                text(-900,-200,['p = ' num2str(stats(2,3))])
                % Exp, t1/t2
                subplot(2,2,3)
                [~,~,stats(3,:)] = getQQdata(x1,x2);
                xlim(alims),ylim(alims)
                title('X1-X2','FontWeight','bold','FontSize',12)
                text(-900,-100,['R^2 = ' num2str(stats(3,1))])
                text(-900,-200,['p = ' num2str(stats(3,3))])
                % Ins, t1/t2
                subplot(2,2,4)
                [~,~,stats(4,:)] = getQQdata(y1,y2);
                xlim(alims),ylim(alims)
                title('Y1-Y2','FontWeight','bold','FontSize',12)
                text(-900,-100,['R^2 = ' num2str(stats(4,1))])
                text(-900,-200,['p = ' num2str(stats(4,3))])
            end
        else
            error('Invalid MAT file loaded!')
        end
    else
        error('Need exactly two files!')
    end
end

function [m,b,stats] = getQQdata(x,y)

% Create quantile-quantile plot
h = qqplot(x,y);

% Get quantiled data
tx = get(h(1),'XData');
ty = get(h(1),'YData');

% Get line
lx = get(h(3),'XData');
ly = get(h(3),'YData');
m = diff(ly)/diff(lx);
b = ly(1) - m*lx(1);
% tl = m*tx + b;
tl = tx;
hold on
plot([tx(1);tx(end)],[tl(1);tl(end)],'-k')
hold off

% Calculate stats
[~,~,~,~,stats] = regress(ty',[tl',ones(length(tl),1)]);


