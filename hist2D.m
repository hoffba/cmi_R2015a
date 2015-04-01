function ha = hist2D(vals,vmin,vmax,nbins)
% Inputs:
%   vals = [N x 2] containing values to histogram
%   vmin = minimum value to bin
%   vmax = maximum value to bin
%   nbins = number of bins

if length(vmin)==1
    vmin = vmin*[1,1];
end
if length(vmax)==1
    vmax = vmax*[1,1];
end
if length(nbins)==1
    nbins = nbins*[1,1];
end

dbin = (vmax-vmin)./nbins;
vmin = vmin + dbin/2;
vmax = vmax - dbin/2;
H = hist3(vals,{linspace(vmin(1),vmax(1),nbins(1)),...
                linspace(vmin(2),vmax(2),nbins(2))});
figure,hi = imagesc(H','XData',[vmin(1),vmax(1)],'YData',[vmin(2),vmax(2)]);
ha = get(hi,'Parent');
axis(ha,'image');
set(ha,'Ydir','normal');
set(get(ha,'YLabel'),'String','Inspiration (HU)')
set(get(ha,'XLabel'),'String','Expiration (HU)')