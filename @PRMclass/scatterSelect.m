% PRMclass function
function stat = scatterSelect(self,imask,tvals,nres)
% Outputs: tmask = logical [nres x nres] mask array of selection
%          ctrs = {1x2} cell array of mask grid locations for x and y

stat = false;
if self.check && ishandle(self.hfscatter) && ...
        strcmp(get(self.hfscatter,'Name'),'PRMscatter')
    
    % Determine selection resolution:
    if (nargin==1)
        nres = 200;
    else
        nres = round(nres);
    end
    
    % Create phantom image in scatterplot axes for masking:
    hold(self.hascatter,'on');
    timg = false(nres);
    xlim = get(self.hascatter,'XLim');
    ylim = get(self.hascatter,'YLim');
    padx = diff(xlim)/(2*nres);
    pady = diff(ylim)/(2*nres);
    hi = image(timg,'Parent',self.hascatter,'AlphaData',timg,...
                    'Xdata',xlim,'Ydata',ylim);
    hold(self.hascatter,'off');
    
    % Create mask from user ROI:
    hfh = imfreehand(self.hascatter);
    tmask = hfh.createMask;
    delete(hi); delete(hfh);
    
    if ~isempty(tmask)
        % Determine image voxel bins:
        [~,bin1] = histc(tvals(:,1),linspace(xlim(1)+padx,xlim(2)-padx,nres));
        [~,bin2] = histc(tvals(:,2),linspace(ylim(1)+pady,ylim(2)-pady,nres));

        % Find selected voxels:
        mskind = find(imask);
        [sy,sx] = ind2sub([nres,nres],find(tmask));
        sel = ismember([bin1,bin2],[sx,sy],'rows');
        nmask = false(self.dims(1:3));
        nmask(mskind(sel)) = true;

        % Replace mask with selection:
        self.mask.merge('replace',nmask);
        
        stat = true;
    end
end