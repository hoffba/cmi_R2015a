% CMIclass function
function genStats(self,~,~)
% Calculates the following VOI statistics:
%       Mean Value
%       Standard Deviation
%       Median Value
%       Volume
%       # Voxels

% only perform stats on loaded image using current mask
if self.img.check && self.img.mask.check 
    [vMean,vStD,vMed,vVol,nVox, VOIstats] = self.img.getStats(self.vec);
    if length(vMean) == 1
    msgbox({['Sum = ',num2str(vMean.*nVox)],...
            ['Mean = ',num2str(vMean)],...
            ['StDev = ',num2str(vStD)],...
            ['Median = ',num2str(vMed)],...
            ['#Voxels = ',num2str(nVox)],...
            ['Volume = ',num2str(vVol)]},...
                    'VOI Stats:');
    else
        fig1 = uifigure('Name','VOI Stats: Editable Table in Workspace');
        fig1.Position = [500 500 370 370];
        uit =uitable(fig1,'Data',VOIstats);
    end
end