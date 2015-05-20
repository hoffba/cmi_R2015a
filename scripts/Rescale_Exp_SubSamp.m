% CMI script
function Rescale_Exp_SubSamp(RegObj)
%   Input: cmiObj = CMIclass object containing current settings
 a=RegObj.cmiObj(1).img.voxsz;a(3)=a(3)*2;RegObj.cmiObj(1).img.setVoxSz(a);
