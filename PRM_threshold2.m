function [slope YInt CI95] = PRM_threshold2(Pixels)
%
% Pixels=cat(1,Pixels1,[Pixels1(:,2) Pixels1(:,1)]);
% size(Pixels)
[coeff,score,latent]=princomp(Pixels);
% coeff is the eigenvector
% score is the new data on new axes
% latent is the eigenvalues.

mx = 5000;
np = size(Pixels,1);
subind = 1:ceil(np/mx):np;

h=figure(200);
scrsz = get(groot,'ScreenSize');
set(h,'Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
subplot(2,2,1)
[p,ErrorEst] = polyfit(Pixels(:,1),Pixels(:,2),1);
[pop_fit,delta] = polyval(p,Pixels(:,1),ErrorEst);
CI=2*delta;
cdate=Pixels(subind,1);pop=Pixels(subind,2);popf=pop_fit(subind);
dscatter(cdate,pop);hold on
plot(cdate,popf,'g-',...
     cdate,popf+CI(subind),'r:',...
     cdate,popf-CI(subind),'r:');hold off
% plot(cdate,pop,'.',...
%      cdate,popf,'g-',...
%      cdate,popf+CI(subind),'r:',...
%      cdate,popf-CI(subind),'r:');
title(['Orig-Data: Slope: ',num2str(p(1,1)),'  Y-Int: ',num2str(p(1,2))...
    ,'  95%CI: ',num2str(mean(CI))]);
xlabel('Pre-Tx Parameter');
ylabel('Post-Tx Parameter');
axis([min(Pixels(:,1)),max(Pixels(:,1)),min(Pixels(:,2)),max(Pixels(:,2))])

subplot(2,2,2)
[p,ErrorEst] = polyfit(score(:,1),score(:,2),1);
[pop_fit,delta] = polyval(p,score(:,1),ErrorEst);
CI=2*delta;
cdate=score(:,1);pop=score(subind,2);
dscatter(cdate(subind),pop);hold on;
plot(cdate(subind),pop_fit(subind),'g-',...
     cdate(subind),pop_fit(subind)+CI(subind),'r:',...
     cdate(subind),pop_fit(subind)-CI(subind),'r:');hold off;
% plot(cdate(subind),pop,'.',...
%      cdate(subind),pop_fit(subind),'g-',...
%      cdate(subind),pop_fit(subind)+CI(subind),'r:',...
%      cdate(subind),pop_fit(subind)-CI(subind),'r:');
title(['Trans-Data: Slope: ',num2str(p(1,1)),'  Y-Int: ',num2str(p(1,2))...
    ,'  95%CI: ',num2str(mean(CI))]);
xlabel('Principal Eigenvector');
ylabel('Second Eigenvector');

Pixels_orig=permute(((coeff)*permute(score,[2,1]))+permute(ones(np,1)*mean(Pixels,1),[2,1]),[2,1]); % This converts data back to original.
Pixels_CIT=permute(((coeff)*permute(cat(2,cdate,pop_fit+CI),[2,1]))+permute(ones(np,1)*mean(Pixels,1),[2,1]),[2,1]); %convert Top CI
Pixels_CIB=permute(((coeff)*permute(cat(2,cdate,pop_fit-CI),[2,1]))+permute(ones(np,1)*mean(Pixels,1),[2,1]),[2,1]); % convert Bottom CI
Pixels_line=permute(((coeff)*permute(cat(2,cdate,pop_fit),[2,1]))+permute(ones(np,1)*mean(Pixels,1),[2,1]),[2,1]); % convert Line

[p_line,ErrorEst] = polyfit(Pixels_line(:,1),Pixels_line(:,2),1); % calculate slope and y-intercept from converted Line
[p_CIT,ErrorEst] = polyfit(Pixels_CIT(:,1),Pixels_CIT(:,2),1); % calculate slope and y-intercept from converted Line
[p_CIB,ErrorEst] = polyfit(Pixels_CIB(:,1),Pixels_CIB(:,2),1); % calculate slope and y-intercept from converted Line

subplot(2,2,3)
dscatter(Pixels_orig(subind,1),Pixels_orig(subind,2));hold on
plot(Pixels_line(subind,1),Pixels_line(subind,2),'r-'...
    ,Pixels_CIT(subind,1),Pixels_CIT(subind,2),'k-'...
    ,Pixels_CIB(subind,1),Pixels_CIB(subind,2),'k-');hold off;
% plot(Pixels_orig(subind,1),Pixels_orig(subind,2),'.',Pixels_line(subind,1),Pixels_line(subind,2),'r-'...
%     ,Pixels_CIT(subind,1),Pixels_CIT(subind,2),'k-'...
%     ,Pixels_CIB(subind,1),Pixels_CIB(subind,2),'k-');
title(['Orig-Data (from Trans): Slope: ',num2str(p_line(1,1)),'  Y-Int: ',num2str(p_line(1,2))...
    ,'  95%CI: ',num2str((p_CIT(1,2)-p_CIB(1,2))/2)]);
% added by MB
slope=p_line(1,1); YInt=p_line(1,2); CI95=(p_CIT(1,2)-p_CIB(1,2))/2;
%%%% added by BL
CI_value = num2str((p_CIT(1,2)-p_CIB(1,2))/2);
assignin('caller','Individual_CI', CI_value);
%%%%

xlabel('Pre-Tx Parameter');
ylabel('Post-Tx Parameter');
axis([min(Pixels(:,1)),max(Pixels(:,1)),min(Pixels(:,2)),max(Pixels(:,2))])