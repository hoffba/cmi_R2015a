% ImageClass function
function stat = autoScale(self,vec)
% Automatically scales CT images based on histogram
stat = false;
if self.check
    
        airHU = -990;
        stHU = 40;
        
    timg = self.mat(:,:,:,vec);
    mask = ~imclose(imclearborder(timg < 0.5*mean(timg(timg>min(timg(:)))),...
                                  cat(3,zeros(3),ones(3),zeros(3))),...
                    strel('disk',round(size(timg,1)/100)));
                
    self.mask.merge('Replace',~mask);
                
    n = timg(mask);
    n(n==min(n)) = []; % removes values outside recon FOV (for clinical data)
    tmin = min(n(:));
    tmax = max(n(:));
    delta = tmax - tmin;
    
    hf = figure('Position',[200,400,1000,400]);
    bins = linspace(tmin,tmax,100);
    h = histc(n(:),bins);
    subplot(1,3,1),plot(bins,h);
    
    % Find Air peak
    bins = linspace(tmin,tmin+0.1*delta,200);%tmin:(tmin+200);
    h = smooth(medfilt1(histc(n(:),bins),20),10);
    h(1) = 0;
    [pks,locs] = findpeaks(h);
    airAU = bins(locs(find(pks==max(pks),1)));
    figure(hf),subplot(1,3,2),plot(bins,h),title(['AUair = ',num2str(airAU)]);
    
    % Find Soft-Tissue peak
    bins = linspace(tmin+0.2*delta,tmin+0.5*delta,200);%(tmin+750):(tmin+1300);
    h = smooth(medfilt1(histc(n(:),bins),20),10);
    [pks,locs] = findpeaks(h);
    stAU = bins(locs(find(pks==max(pks),1)));
    figure(hf),subplot(1,3,3),plot(bins,h),title(['AUst = ',num2str(stAU)]);
    
    m = (stHU - airHU) / (stAU - airAU);
    b = -m * airAU + airHU;
    
    self.imgScale( vec ,...
                  m * self.scaleM(vec) ,...
                  m * self.scaleB(vec) + b);
    stat = true;
end