% ImageClass function
% Generates a P-P plot and returns statistics D and AUC
%   - Must have at least two images and a mask loaded
%   - Analyzes first two images (#1=X and #2=Y)
%   - Mask is thresholded using current image thresholds
% Input:  opts = (1) Wiener filter : 1=yes, 0=no
%                (2) PP Plot : 1=yes, 0=no
%                (3) X-Histogram Fit : 1=all, 0=no
%                (4) Y-Histogram Fit : 2=Max, 1=all, 0=no
%                (5) PCA :  1=yes, 0=no
% Output: stats - structure containing analysis statistics
%                   (also displayed in popup window)
function stats = ppPlot(self,cvec,opts)

stats = [];
if self.check && (self.dims(4)>1) && self.mask.check && (nargin>1)
    dvec = self.prm.dvec; dvec(dvec==0) = cvec;
    
    try
        lims = self.prm.cutoff; lims(lims(:,1)==0) = cvec;
    catch
        lims = cat(1,[1 2],self.thresh')'; lims(lims(:,1)==0) = cvec;
    end
    
    xi = find(lims(:,1)==dvec(1),1);
    yi = find(lims(:,1)==dvec(2),1);
    
    if (nargin<3) || ~isnumeric(opts) || (length(opts)~=5)
        % Default options:
        opts = [ 1 0 0 2 1 ];
    end
    
    if ~(isempty(xi) || isempty(yi)) && (xi~=yi)
        
        hw = waitbar(0,'Grabbing thresholded mask and image values ...');
        tmask = self.getThreshMask(dvec);
        doff = prod(self.dims(1:3));
        ii = find(tmask);
        
        % Wiener filter for de-noising:
        if opts(1)==1
            self.imgFilt(dvec,'median',[3,3]);
        end

        % Extract image values for analysis:
        xvals = self.mat(ii+doff*(dvec(1)-1));
        yvals = self.mat(ii+doff*(dvec(2)-1));

        % PP-plot analysis:
        if opts(2)==1
            if (~isfield(self.h,'ppfig') || ~ishandle(self.h.ppfig))
                self.h.ppfig = figure('Name','PPplot'); self.h.ppax = axes;
            end
            waitbar(1/4,hw,'Generating P-P plot and statistics ...');
            [~,stats.D,stats.AUC] = PPplot(xvals,yvals,self.h.ppax);
        end
        
        % X-histogram analysis:
        if opts(3)==1
            waitbar(2/4,hw,'Determining X-histogram confidence intervals ...');
            if ~isfield(self.h,'histfitX') || ~ishandle(self.h.histfitX)
                self.h.histfitX = figure('Name','Histogram Fit (X)');
            else
                figure(self.h.histfitX);
            end
            [stats.xEsts,stats.xHist,stats.xCI,stats.xGoF] = histo_fit(xvals,lims(xi,2),lims(xi,3),4); % 4 = trunc GEV
        end
        
        % Y-histogram analysis:
        if ismember(opts(4),1:2)
            waitbar(3/4,hw,'Determining Y-histogram confidence intervals ...');
            if ~isfield(self.h,'histfitY') || ~ishandle(self.h.histfitY)
                self.h.histfitY = figure('Name','Histogram Fit (Y)');
            else
                figure(self.h.histfitY);
            end
            % CJG added to calculate histogram through max of JDH
            % 29Oct2015 (based on JDH_bounds_v3.m)
            if opts(4)==2
                nbins = 50;
                [N,C] = hist3([xvals yvals],[nbins nbins]);
                ma = max(max(N));
                [ii,~] = ind2sub(size(N),find(ma==N)); % This identifies were the max value of the JDH is located in [i,j]
                ii = ii(1,1);              % J=J1(1,1);
                h1 = (C{1}(2)-C{1}(1))/2; % h2=(C{2}(2)-C{2}(1))/2;
                %             exp_temp=expv(insv<(C{2}(J)+h2)&insv>=(C{2}(J)-h2)); % use the pre-filtered data insv and expv
                ind = xvals<(C{1}(ii)+h1) & xvals>=(C{1}(ii)-h1);
                disp(['Max of JDH: ',num2str(ii)])
            else
                ind = true(size(yvals));
            end
            [stats.yEsts,stats.yHist,stats.yCI,stats.yGoF] = histo_fit(yvals(ind),lims(yi,2),lims(yi,3),4);
        end
        
        % PCA analysis:
        if opts(5)==1
            waitbar(4/4,hw,'Determining Principal Component Analysis ...');
            [stats.PCA.coeff,stats.PCA.score,stats.PCA.latent] = pca([xvals yvals]);
            stats.PCA.theta = round(45-(atand(stats.PCA.coeff(1,1)/stats.PCA.coeff(2,1))),2);
        end
        delete(hw);
        if ~nargout
            assignin('base','PPstats',stats);
        end
    end
end
