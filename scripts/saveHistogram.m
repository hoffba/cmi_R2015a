% CMI script
function saveHistogram(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings


 
 %% plot histogram with bars
%  og = [1 0.5 0];
%  bl = [0 0 1];
%  cmiObj.setVec(1);
%  AbinVals=cmiObj.histObj.binVals;
%  AbinLocs=cmiObj.histObj.binLocs;
%  figure(50);
%  bar(AbinLocs,AbinVals(:,1),'FaceColor',og,'EdgeColor',og); % Full histogram
%  
%  cmiObj.setVec(2);
%  AbinVals=cmiObj.histObj.binVals;
%  AbinLocs=cmiObj.histObj.binLocs;
%  figure(50);
%  hold on;
%  bar(AbinLocs,AbinVals(:,1),'FaceColor',bl,'EdgeColor',bl); % Slice
%  hold off;

 %% plot histogram with lines
 og = [1 0.5 0];
 bl = [0 0.3 1];
 gd = [0.6 0.8 0];
 p = [1 0 1];
 g = [0 1 0];
 mnt = [0 1 0.7];
 
 cs = [og; bl; gd; mnt; p; g];
 
 cmiObj.setHistVis(1);
 
 for i=1:cmiObj.img.dims(4)
     cmiObj.setVec(i);
     AbinVals=cmiObj.histObj.binVals;
     AbinLocs=cmiObj.histObj.binLocs;
%      hmean(i) = mean(AbinVals(:,1));
%      hmed(i) = median(AbinVals(:,1));
     figure(50);
     plot(AbinLocs,AbinVals(:,1),'Color',cs(i,:)); % Full histogram 1
     hold on;
 end
 
hold off;

%% histogram stats
% fprintf('HISTOGRAM STATS for %d vectors\n',cmiObj.img.dims(4));
% fprintf('===================================\n');
% fprintf('Hist mean=%d\n',hmean);
% fprintf('Hist median=%d\n',hmed);

end

 