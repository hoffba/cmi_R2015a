% ImageClass function
function [CI,slp] = calcCI(self,vec)

if self.check && self.mask.check && (self.dims(4)>1) && (length(vec)==2)
    
    vals = self.getMaskVals(vec);
    vmeans = mean(vals,1);
    ext = [ min(vals,[],1) ; max(vals,[],1) ];
    np = size(vals,1);
    maxplotn = 5000;
    npsub = min(np,maxplotn);
    intvl = ceil(np/maxplotn);
    cmap = jet(64); cmap(1,:) = [1,1,1];
    xlocs = linspace(ext(1,1),ext(2,1),128);
    ylocs = linspace(ext(1,2),ext(2,2),128);
    
    figure;
    subplot(2,2,1)
    [p,ErrorEst] = polyfit(vals(:,1),vals(:,2),1);
    [pop_fit,delta] = polyval(p,ext(:,1),ErrorEst);
    CI=2*mean(delta);
    H = hist2(vals(:,1),vals(:,2),xlocs,ylocs);
    figure,imagesc(H); colormap(cmap);
    set(gca,'YDir','normal');
    plotlines(gca,p(1)*[1,1,1],[p(2),p(2)+CI,p(2)-CI]/diff(ext(:,1))*128);
    hold on,
    plot( ext(:,1) , pop_fit , 'g-',...
          ext(:,1) , pop_fit+2*delta , 'r:',...
          ext(:,1) , pop_fit-2*delta , 'r:');
    title(['Orig-Data: Slope: ',num2str(p(1,1)),'  Y-Int: ',num2str(p(1,2))...
        ,'  95%CI: ',num2str(CI)]);
    xlabel('Pre-Tx Parameter');
    ylabel('Post-Tx Parameter');
    axis(ext(:)')

    [coeff,vals,~] = princomp(vals);
    % % coeff is the eigenvector
    % % score is the new data on new axes
    % % latent is the eigenvalues.
    subplot(2,2,2)
    [p,ErrorEst] = polyfit(vals(:,1),vals(:,2),1);
    [pop_fit,delta] = polyval(p,ext(:,1),ErrorEst);
    CI = 2*mean(delta);
    %cdate=score(:,1);pop=score(:,2);
    plot( vals(1:intvl:end,1) , vals(1:intvl:end,2) , '.',...
          ext(:,1) , pop_fit , 'g-',...
          ext(:,1) , pop_fit+2*delta , 'r:',...
          ext(:,1) , pop_fit-2*delta , 'r:');
    title(['Trans-Data: Slope: ',num2str(p(1,1)),'  Y-Int: ',num2str(p(1,2))...
        ,'  95%CI: ',num2str(CI)]);
    xlabel('Principal Eigenvector');
    ylabel('Second Eigenvector');
    
%     np = size(vals,1); % sub-sampled # of values
    Pixels_orig = (vals(1:intvl:end,:)*coeff'    + (ones(npsub,1)*vmeans) ); % This converts data back to original.
    Pixels_CIT =  ([ext(:,1),pop_fit+CI]*coeff'  + ([1,1]*vmeans) ); % convert Top CI
    Pixels_CIB =  ([ext(:,1),pop_fit-CI]*coeff'  + ([1,1]*vmeans) ); % convert Bottom CI
    Pixels_line = ([ext(:,1),pop_fit]*coeff'     + ([1,1]*vmeans) ); % convert Line
    
    [p_line,~] = polyfit(Pixels_line(:,1),Pixels_line(:,2),1); % calculate slope and y-intercept from converted Line
    [p_CIT,~] = polyfit(Pixels_CIT(:,1),Pixels_CIT(:,2),1); % calculate slope and y-intercept from converted Line
    [p_CIB,~] = polyfit(Pixels_CIB(:,1),Pixels_CIB(:,2),1); % calculate slope and y-intercept from converted Line
    CI = (p_CIT(1,2)-p_CIB(1,2))/2;
    slp = p_line(1,1);
    subplot(2,2,3)
    plot( Pixels_orig(:,1) , Pixels_orig(:,2) , '.',...
          Pixels_line(:,1) , Pixels_line(:,2) , 'r-',...
          Pixels_CIT(:,1)  , Pixels_CIT(:,2)  , 'k-',...
          Pixels_CIB(:,1)  , Pixels_CIB(:,2)  , 'k-');
    title(['Orig-Data (from Trans): Slope: ',num2str(slp),'  Y-Int: ',num2str(p_line(1,2))...
        ,'  95%CI: ',num2str(CI)]);
    xlabel('Pre-Tx Parameter');
    ylabel('Post-Tx Parameter');
    axis(ext(:)')
end