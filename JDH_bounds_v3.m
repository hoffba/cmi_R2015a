function PRMv3=JDH_bounds_v3(data);


if nargin==0
    cmiObj=evalin('base','cmiObj0')
    exp=cmiObj.img.mat(:,:,:,1);ins=cmiObj.img.mat(:,:,:,2);
    mask=cmiObj.img.mask.mat;
    expv=exp(mask);insv=ins(mask);
else
    exp=data{1};ins=data{2};mask=data{3};
    expv=exp(mask);insv=ins(mask);
end
VOI=numel(mask);

UL=-250;
exp1=expv(expv<=UL&insv<=UL);
ins1=insv(expv<=UL&insv<=UL);
assignin('base','CT_exp_ins',[exp1, ins1]);
%% PPplot
'PPplot'
figure(101); 
[~,stats.D,stats.AUC] = PPplot(exp1,ins1,gca);
%%
% Principal Component Analysis
'Principal Component Analysis'
[coeff,score,latent]=pca([exp1 ins1]);
theta=45-(atand(coeff(2,1)/coeff(1,1)))

%% 3D Histogram
'3D Histogram'
func_tag=4;
k=2;

nbins=50;
CI=zeros(2,9);
goF=zeros(1,1);

[N,C]=hist3([exp1 ins1],[nbins nbins]);

figure(99);h_hist=imagesc(C{1},C{2},N);
hold on
plot(-1024:100:0,-1024:100:0,'-w')
hold off
gca;axis([-1000 -300 -1000 -300]);
set(gca,'YDir','normal');

assignin('base','JDH',[{C{1}},{C{2}},{N}]);
% figure(100);imagesc(Ca{1},Ca{2},Na);

% [a1,a2,a3]=pca(N);
a2=N;
ma=max(max(a2));
[I1,J1]=ind2sub(size(a2),find(ma==a2)); % This identifies were the max value of the JDH is located in [i,j]
I=I1(1,1);
J=J1(1,1);
assignin('base','CenterJDH',[I,J]);

%% Generate histogram from bins
% countI=1;countJ=1;
% for i=1:size(N,2)
%     if ~isempty(ones(N(I,i),1)) % Only iterate through row I ( you don't use the whole histogram just the histogram through the peak value)
%         if countI==1
%             mY=(C{1}(i)*ones(N(I,i),1)); % Inspiration
%         else
%             mY=cat(1,mY,C{1}(i)*ones(N(I,i),1));
%         end
%         climI(countI)=C{1}(i);countI=countI+1;
%     end
%      if ~isempty(ones(N(i,J),1)) % Only iterate through column J
%         if countJ==1
%             mX=(C{2}(i)*ones(N(i,J),1)); % Excitation
%         else
%             mX=cat(1,mX,C{2}(i)*ones(N(i,J),1));
%         end
%         climJ(countJ)=C{2}(i);countJ=countJ+1;
%     end
% end

% grabbed code from cjg_fitfunc_v5_proc and histo_fit
% for q=1:2
%     if q==1
%         'GEV'
%         %         mXY=mY;clim=climI;
%         mXY=ins_temp;
%         
%     else
%         %         mXY=mX;clim=climJ;
%         mXY=exp_temp;
%     end
%     
%     figure(10+q);
%     [paramEsts,paramCIs,goodF]=histo_fit(mXY,min(mXY),max(mXY),func_tag);
%     %         [paramEsts,paramCIs,goodF]=histo_fit(mXY,min(clim),max(clim),func_tag);
%     goodF
%     CI(q,1:9)=round(paramEsts(4:end));goF(q,1)=goodF;
%     GEV_quant(q,1:12)=round(paramEsts(1:end));
% end
% 
%% GEV for Exp(J) and Ins(I) histograms
'GEV for Exp(J) and Ins(I) histograms'
h1=(C{1}(2)-C{1}(1))/2;h2=(C{2}(2)-C{2}(1))/2;
exp_temp=expv(insv<(C{2}(J)+h2)&insv>=(C{2}(J)-h2)); % use the pre-filtered data insv and expv
ins_temp=insv(expv<(C{1}(I)+h1)&expv>=(C{1}(I)-h1));

% Inspriation
figure(11)
[paramEsts,paramCIs,goodF]=histo_fit(ins_temp,-1000,-250,func_tag);
CI(1,1:9)=round(paramEsts(4:end),4);goF(1,1)=goodF;
GEV_quant(1,1:12)=round(paramEsts(1:end),4);
% Expiration
figure(12)
[paramEsts,paramCIs,goodF]=histo_fit(exp_temp,-1000,-250,func_tag);
CI(2,1:9)=round(paramEsts(4:end),4);goF(2,1)=goodF
GEV_quant(2,1:12)=round(paramEsts(1:end),4);
%%
'GEV for Exp(all) and Ins(all) histograms'
% Inspriation
figure(13)
[paramEsts,paramCIs,goodF]=histo_fit(ins1(1:1000:end),-1000,-250,func_tag);
GEV_quant2(1,1:12)=round(paramEsts(1:end),4);
GoF2(1,1)=goodF;
% Expiration
figure(14)
[paramEsts,paramCIs,goodF]=histo_fit(exp1(1:1000:end),-1000,-250,func_tag);
GEV_quant2(2,1:12)=round(paramEsts(1:end),4);
GoF2(2,1)=goodF
%% This is to determine the JDH values for each quantile
% quantiles are 0.025 0.05 0.15 0.25 0.5. 
nN=N/max(max(N));
[X,Y]=meshgrid(C{2},C{1});
fnN=6*ones(size(nN));

for i=5:-1:1
%     J(i,1)=ind2sub(size(C{1}),find(CI(1,i)>=C{1}));
    nJ(i,1)=max(ind2sub(size(C{2}),find(CI(1,i)>=C{2})));
    nNI(i,1)=nN(I,nJ(i,1));
    fnN(nN<=nNI(i,1))=i;
end
% This figure generates the JDH using the quantils
% figure(102);imagesc(C{1},C{2},fnN);
%%
% Save results to workspace
% PRMv3=[VOI,Upper Limit,D,AUC,Latent,Coeff,theta,Exp_JDHmax,Ins_JDHmax,paramEsts_Ins,Gof_Ins,paramEsts_Exp,Gof_Exp,...
% paramEsts_Ins2,Gof_Ins2,paramEsts_Exp2,Gof_Exp2]
PRMv3=cat(1,[VOI,UL,stats.D,stats.AUC,reshape(latent,[1,2]),reshape(coeff,[1,4]),theta,...
    C{2}(I),C{1}(J),GEV_quant(1,1:end),goF(1,1),GEV_quant(2,1:end),goF(2,1),mean(exp1),mean(ins1),std(exp1),std(ins1),...
    GEV_quant2(1,1:end),GoF2(1,1),GEV_quant2(2,1:end),GoF2(2,1)]);
assignin('base','PRMv3',PRMv3);

assignin('base','GEV_results_Ins',cat(1,[goF(1,1),CI(1,1:5),nNI']));

