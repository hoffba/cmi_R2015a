function [stats,D,AUC] = PPplot(x,y,h)
% Generates a Probability-Probability plot
% Returns statistics comparing the generated plot to unity
if (nargin>1)
    [Fx,xx] = ecdf(x(:));
    [Fy,xy] = ecdf(y(:));
    Fx(1)=[];xx(1)=[];Fy(1)=[];xy(1)=[];

    % x-axis must be aligned
    xmin = min(xx(1),xy(1));
    xmax = max(xx(end),xy(end));
    xi = (xmin:xmax)';
    Fxi = interp1(xx,Fx,xi);
    Fyi = interp1(xy,Fy,xi);
    n = length(xi);

    % Generate Statistics
    [~,~,~,~,stats] = regress(Fyi,[Fxi,ones(n,1)]);
    % stats = [ R^2 , F , p-value , error(variance)]
    D = max(abs(Fyi-Fxi)) + 1/(2*n);
    Fxii=Fxi;Fyii=Fyi;
    Fxii(isnan(Fxi)|isnan(Fyi))=[];Fyii(isnan(Fxi)|isnan(Fyi))=[];
    AUC=simps(Fxii,Fyii,1)-0.5;
    %AUC(isnan(AUC))=1000;
    % Dks = max(abs(Fyi - Fxi));
    % alpha = 0.05;
    % K = kolmcdf(sqrt((nx+ny)/(nx*ny))*abs(Fxi-Fyi));
    % Ka = xx(abs(K-1+alpha)==min(abs(K-1+alpha)));
    % h = sqrt(nx*ny/(nx+ny))*Dks > Ka;

    % Create P-P Plot
    if (nargin<3) || isempty(h) || ~ishandle(h)
        figure, h = axes;
    end
    plot(h,Fxi,Fyi,'.b',[0,1],[0,1],'--r',[0,1-D],[D,1],'g',[D,1],[0,1-D],'g')
    %text(0.15,0.85,['R^2 = ' num2str(stats(1))],'FontSize',12)
    text(0.7,0.2,['D = ' num2str(D)],'FontSize',9,'FontAngle','italic','Parent',h) 
    text(0.7,0.1,['AUC = ' num2str(AUC)],'FontSize',9,'FontAngle','italic','Parent',h)
    %title(h,'P-P plot','FontSize',12,'FontWeight','bold')
    xlabel(h,'Fx','FontSize',12),ylabel(h,'Fy','FontSize',12)
end

