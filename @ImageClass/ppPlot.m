% ImageClass function
% Generates a P-P plot and returns statistics D and AUC
%   - Must have at least two images and a mask loaded
%   - Analyzes first two images (#1=X and #2=Y)
%   - Mask is thresholded using current image thresholds
% Output: stats - structure containing analysis statistics
%                   (also displayed in popup window)
function stats = ppPlot(self,cvec)

stats = [];
if self.check && (self.dims(4)>1) && self.mask.check
    % Get PRM settings:
    if isa(cvec,'CMIclass')
        cvec = cvec.vec;
    end
    dvec = self.prm.dvec; dvec(dvec==0) = cvec;
    lims = self.prm.cutoff; lims(lims(:,1)==0) = cvec;
    xi = find(lims(:,1)==dvec(1),1);
    yi = find(lims(:,1)==dvec(2),1);
    
    if ~(isempty(xi) || isempty(yi)) && (xi~=yi)
        hw = waitbar(0,'Grabbing thresholded mask and image values ...');
        tmask = self.getThreshMask(dvec(1)) & self.getThreshMask(dvec(2));
        doff = prod(self.dims(1:3));
        ii = find(tmask);
        xvals = self.mat(ii+doff*(dvec(1)-1));
        yvals = self.mat(ii+doff*(dvec(2)-1));
        if ~isfield(self.h,'ppfig') || ~ishandle(self.h.ppfig)
            self.h.ppfig = figure('Name','PPplot'); self.h.ppax = axes;
        end
        waitbar(1/4,hw,'Generating P-P plot and statistics ...');
        [~,stats.D,stats.AUC] = PPplot(xvals,yvals,self.h.ppax);
        waitbar(2/4,hw,'Determining X-histogram confidence intervals ...');
        if ~isfield(self.h,'histfitX') || ~ishandle(self.h.histfitX)
            self.h.histfitX = figure('Name','Histogram Fit (X)');
        else
            figure(self.h.histfitX);
        end
        [stats.xEsts,stats.xCI,stats.xGoF] = histo_fit(xvals,lims(xi,2),lims(xi,3),4); % 4 = trunc GEV
        waitbar(3/4,hw,'Determining Y-histogram confidence intervals ...');
        if ~isfield(self.h,'histfitY') || ~ishandle(self.h.histfitY)
            self.h.histfitY = figure('Name','Histogram Fit (Y)');
        else
            figure(self.h.histfitY);
        end
        [stats.yEsts,stats.yCI,stats.yGoF] = histo_fit(yvals,lims(yi,2),lims(yi,3),4);
        delete(hw);
        if ~nargout
            assignin('base','PPstats',stats);
        end
    end
end
  