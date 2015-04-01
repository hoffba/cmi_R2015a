% CMI script
function example(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings

disp('Hello World!')

if nargin && isa(cmiObj,'CMIclass')
    disp(['Image size: ',num2str(cmiObj.img.dims)])
end