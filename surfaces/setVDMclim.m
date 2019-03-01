function setVDMclim(hf,r)

ha = findobj(hf.Children,'Type','Axes');
hc = findobj(hf.Children,'Type','Colorbar');

tics = linspace(-r,r,6);
ticlab = num2str(exp(tics)',2);
set(hc,'Limits',[-r,r],'Ticks',tics,'TickLabels',ticlab);

set(ha,'CLim',[-r,r]);