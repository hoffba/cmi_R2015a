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
    if isempty(lims)
        lims = cat(1,[1 2],self.thresh')'; lims(lims(:,1)==0) = cvec;
        if dvec(1)==dvec(2)
            dvec(2)=dvec(1)+1;
        end
    end
    
    xi = find(lims(:,1)==dvec(1),1);
    yi = find(lims(:,1)==dvec(2),1);
    
    prompt={'Wiener2[3 3]: 1=yes, 0=no','PP Plot: 1=yes, 0=no','X-Histogram Fit: 1=all, 0=no','Y-Histogram Fit: 2=Max, 1=all, 0=no',...
        'PCA: 1=yes, 0=no','Stanford Threshold: 1=yes, 0=no'};
    dlg_title='Which Analysis?';
    num_lines=1;
    def={'0','0','0','0','0','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if ~(isempty(xi) || isempty(yi) || isempty(answer)) && (xi~=yi)
        
        hw = waitbar(0,'Grabbing thresholded mask and image values ...');
        tmask = self.getThreshMask(dvec(1)) & self.getThreshMask(dvec(2));
        doff = prod(self.dims(1:3));
        ii = find(tmask);
        
        
        ximg=self.mat(:,:,:,1);yimg=self.mat(:,:,:,2);
        if strcmp(answer{1},'1')
            hing = waitbar(0,'Wiener2 Filter ...');
            for i=1:size(ximg,3)
                waitbar(i/size(ximg,3),hing)
                ximg(:,:,i)=wiener2(ximg(:,:,i),[3 3]);
                yimg(:,:,i)=wiener2(yimg(:,:,i),[3 3]);
            end
            close(hing);
        end  
        xvals = ximg(ii+doff*(dvec(1)-1));
        yvals = yimg(ii+doff*(dvec(1)-1));
%         xvals = self.mat(ii+doff*(dvec(1)-1));
%         yvals = self.mat(ii+doff*(dvec(2)-1));

        if strcmp(answer{2},'1')
            if (~isfield(self.h,'ppfig') || ~ishandle(self.h.ppfig))
                self.h.ppfig = figure('Name','PPplot'); self.h.ppax = axes;
            end
            waitbar(1/5,hw,'Generating P-P plot and statistics ...');
            [~,stats.D,stats.AUC] = PPplot(xvals,yvals,self.h.ppax);
        end
        if strcmp(answer{3},'1')
            waitbar(2/5,hw,'Determining X-histogram confidence intervals ...');
            if ~isfield(self.h,'histfitX') || ~ishandle(self.h.histfitX)
                self.h.histfitX = figure('Name','Histogram Fit (X)');
            else
                figure(self.h.histfitX);
            end
            [stats.xEsts,stats.xHist,stats.xCI,stats.xGoF] = histo_fit(xvals,lims(xi,2),lims(xi,3),4); % 4 = trunc GEV
        end
        if strcmp(answer{4},'1')||strcmp(answer{4},'2')
            waitbar(3/5,hw,'Determining Y-histogram confidence intervals ...');
            if ~isfield(self.h,'histfitY') || ~ishandle(self.h.histfitY)
                self.h.histfitY = figure('Name','Histogram Fit (Y)');
            else
                figure(self.h.histfitY);
            end
            % CJG added to calculate histogram through max of JDH
            % 29Oct2015 (based on JDH_bounds_v3.m)
            if strcmp(answer{4},'2')
                nbins=50;
                [N,C]=hist3([xvals yvals],[nbins nbins]);
                a2=N;
                ma=max(max(a2));
                [I1,J1]=ind2sub(size(a2),find(ma==a2)); % This identifies were the max value of the JDH is located in [i,j]
                I=I1(1,1);
                J=J1(1,1);
                h1=(C{1}(2)-C{1}(1))/2;h2=(C{2}(2)-C{2}(1))/2;
                %             exp_temp=expv(insv<(C{2}(J)+h2)&insv>=(C{2}(J)-h2)); % use the pre-filtered data insv and expv
                yvals_temp=yvals(xvals<(C{1}(I)+h1)&xvals>=(C{1}(I)-h1));
                ['Max of JDH: ',num2str(I)]
            else
                yvals_temp=yvals;
            end
            [stats.yEsts,stats.yHist,stats.yCI,stats.yGoF] = histo_fit(yvals_temp,lims(yi,2),lims(yi,3),4);
        end
        if strcmp(answer{5},'1')
            waitbar(4/5,hw,'Determining Principal Component Analysis ...');
            [stats.PCA.coeff,stats.PCA.score,stats.PCA.latent]=pca([xvals yvals]);
            stats.PCA.theta=round(45-(atand(stats.PCA.coeff(2,1)/stats.PCA.coeff(1,1))),2);
        end
        if strcmp(answer{6},'1')
            waitbar(5/5,hw,'Calculate Stanford Thresholds ...');
            % This technique of calculating the gas trapping threshold is
            % published in CHEST/123/5/May,2003/Goris et al.
            % T(i)=X-(i-1)(X-Y)/3-1(1-D/343)(X-Y)/3
            % X = ins 90%
            % Y = ins 50%
            % D = ins90% - exp90%
            x90=round(quantile(xvals,0.9)); % 90%ils from xvals
            y90=round(quantile(yvals,0.9)); % 90%ils from yvals (X in equation)
            D=y90-x90; % change in 90%ile values in ins and exp (D in equation)
            y50=round(quantile(yvals,0.5)); % 50%ils from yvals (Y in equation)
            for i=1:3
                stats.Stanford.T(i)=y90-(i-1)*(y90-y50)/3-(1-D/343)*(y90-y50)/3;
            end
            stats.Stanford.x90=x90;
            stats.Stanford.y90=y90;
            stats.Stanford.D=D;
            stats.Stanford.y50=y50;
        end
        delete(hw);
        if ~nargout
            assignin('base','PPstats',stats);
        end
    end
end
