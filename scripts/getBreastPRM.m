% CMI script
function getBreastPRM(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings

[labels,pcts]=cmiObj.img.prm.getStats;
assignin('base','labels',labels);assignin('base','pcts',pcts)
outvals = [labels; num2cell(pcts)];
xlwrite(fullfile(pwd,'PRMtmp.xls'),outvals);