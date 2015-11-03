% Add lines to a current plot
function plotlines(ha,m,b,linespec)
% Inputs: ha  = handle to axes object
%         m   = vector of line slopes
%         b   = vector of line y-intercepts
if nargin<4
    linespec = repmat({'-k'},1,length(m));
end
if (nargin>=3) && ishandle(ha) && (length(m)==length(b)) && iscellstr(linespec)
    xlim = get(ha,'XLim');
    ylim = get(ha,'YLim');
    hold(ha,'on')
    for i = 1:length(m)
        if isinf(m(i))
            x = b(i)*[1,1];
            y = ylim;
        else
            x = xlim;
            y = m(i)*xlim + b(i);
        end
<<<<<<< HEAD
        plot(ha,x,y,'-k','Tag',['Line',num2str(i)]);
        plot(ha,x,x,'--k','Tag',['Line',num2str(i)]);
=======
        plot(ha,x,y,linespec{i},'Tag',['Line',num2str(i)]);
>>>>>>> origin/master
    end
    hold(ha,'off')
end