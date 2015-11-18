function [paramEsts,Hist_results,paramCIs,goodF]=histo_fit(tdata,LB,UB,func)

% code acquired from
% www.mathworks.com/help/stats/examples/fitting-custom-univariate-distributions.html

if nargin==1||isinf(LB) 
    UB=-200;
    LB=-1000;
    func=4; % 1=bi-Gaussian, 2=trunc Gaussian, 3=trunc Gamma, 4=trunc GEV
end


tdata=tdata(tdata>LB&tdata<UB);
if numel(tdata)>100000
    n=100000;
else
    n=numel(tdata);
end
x=tdata(1:round(numel(tdata)/n):end);
mu=mean(x);
sigma=sqrt(std(x));

% if isempty(x_clim)
    nbins=50;
    h_clim=(UB-LB)/nbins;
    x_clim=LB:h_clim:UB;
% else
%     nbins=numel(x_clim);
% end

[f,xi]=hist(x,x_clim);
% finx=find(f~=0); % This was to remove zeros from the histogram
subplot(1,3,1);plot(xi,f/max(f),'r');axis([LB UB 0 1]);
pause(1)

if func==1
    % normalize CT lung data
    logCT=log(x+1000+1024);
    max_logCT=max(logCT);
    x=logCT/max_logCT;
    [u,sig,t,q]=fit_mix_gaussian(x,2);
    
    % Unnormalize data
    paramEsts(1)=(exp(u.*max_logCT)-1000-1024);
    paramEsts(2)=sig;
    u_plus=u+1.6*sig;u_minus=u-1.6*sig;
    paramEsts(3)=(exp(u_plus.*max_logCT)-1000-1024);
    paramEsts(4)=(exp(u_minus.*max_logCT)-1000-1024);
    paramEsts(5)=0;
elseif func==2
    pdf_truncfunc=@(x,mu,sigma) normpdf(x,mu,sigma)./normcdf(UB,mu,sigma);start=[mu std(x)];
    [paramEsts,paramCIs]=mle(x,'pdf',pdf_truncfunc,'start',start, 'lower', [-Inf 0]);
    y=normpdf(xi,paramEsts(1),paramEsts(2));
    paramEsts(3)=paramEsts(1)+1.6*(paramEsts(2));
    paramEsts(4)=paramEsts(1)-1.6*(paramEsts(2));
    paramEsts(5)=0;
elseif func==3
    x1=x+1024;a=5;b=30;
    pdf_truncfunc=@(x1,a,b) gampdf(x1,a,b)./gamcdf(UB+1024,a,b);start=[a b];
    %     [paramEsts,paramCIs]=gamfit(x+1024);
    [paramEsts,paramCIs]=mle(x1,'distribution','Gamma','pdf',pdf_truncfunc,'start',start, 'lower', [0 Inf]);
    y=gampdf(xi+1024,paramEsts(1),paramEsts(2));
    GR=gamrnd(paramEsts(1),paramEsts(2),[10000,1]);
    Quant_results=quantile(GR(:,1),[0.025 0.975]); % 95% confidence interval [2.5% 97.5%]
    Hist_results=quantile(x,[0.025 0.05 0.15 0.25 0.5 0.75 0.85 0.95 0.975]); % quantiles from data
    paramEsts=cat(2,paramEsts,Quant_results-1024);
    paramEsts(5)=0;
else
    x1=x+1024;a=0;b=100;c=200;
    pdf_truncfunc=@(x1,a,b,c) gevpdf(x1,a,b,c)./gevcdf(UB+1024,a,b,c);start=[a b c];
%     [paramEsts,paramCIs]=gamfit(x+1024);
    [paramEsts,paramCIs]=mle(x1,'distribution','gev','pdf',pdf_truncfunc,'start',start, 'upper', [0 Inf]);
    y=gevpdf(xi+1024,paramEsts(1),paramEsts(2),paramEsts(3));
    GR=gevrnd(paramEsts(1),paramEsts(2),paramEsts(3),[10000,1]);
    Quant_results=round(quantile(GR(:,1),[0.025 0.05 0.1 0.15 0.25 0.5 0.75 0.85 0.9 0.95 0.975])); % quantile from fit
    Hist_results=round(quantile(x,[0.025 0.05 0.1 0.15 0.25 0.5 0.75 0.85 0.9 0.95 0.975])); % quantiles from data
    paramEsts=cat(2,paramEsts,Quant_results-1024);
end

subplot(1,3,2);plot(xi,y/max(y));axis([LB UB 0 1]);
subplot(1,3,3);plot(xi,f/max(f),'r');axis([LB UB 0 1]);hold on;plot(xi,y/max(y));axis([LB UB 0 1]);hold off;

goodF=1-gfit(f/max(f),y/max(y),'4'); %4=NRMSE
% Used to check Gamma fit over gamrand using fit [a b]
% [nelements,centers]=hist(GR,100);
% figure(200);bar(centers-1024,nelements/max(nelements));hold on;
% plot(xi,y/max(y));hold off;





