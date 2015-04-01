% CMI script
function calcJacobian(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings
% useful features for identifying methods and properties
% >>methods ImageClass
% >>properties ImageClass

tic;
'Load X disp field (hit cancel to skip)'
[data(:,:,:,1),label,fov1,fnameOut] = cmi_load([],[],[]);
if ~isempty(data)
    'Load Y disp field'
    [data(:,:,:,2),label,fov1,fnameOut] = cmi_load([],[],[]);
    'Load Z disp field'
    [data(:,:,:,3),label,fov1,fnameOut] = cmi_load([],[],[]);
    'Load Exp_Label'
    [mask(:,:,:),label,fov,fnameOut] = cmi_load([],[],[]);
    
    mask(mask==3)=0;mask(mask>0)=1;mask=logical(mask);
    matrix=size(mask);
    voxsz=fov1./matrix;
    clear fov fnameOut label
else
    'Analyzing Data from CMI GUI'
    data = cmiObj.img.mat(:,:,:,1:3);
    % Normalize voxels
    voxsz = cmiObj.img.voxsz;
    matrix = cmiObj.img.dims(1:3);
    mask = cmiObj.img.mask.mat;
end

% Take gradient of displacement maps
% [dFx,dFy,dFz] = gradient(F,hx,hy,hz);
[uxx,uxy,uxz]=gradient(-data(:,:,:,2)*voxsz(2),voxsz(2),voxsz(1),voxsz(3)); % X displacement map
[uyx,uyy,uyz]=gradient(-data(:,:,:,1)*voxsz(1),voxsz(2),voxsz(1),voxsz(3)); % Y displacement map
[uzx,uzy,uzz]=gradient(-data(:,:,:,3)*voxsz(3),voxsz(2),voxsz(1),voxsz(3)); % Z displacement map

% maxvoxsz=max(voxsz);
% [uxx,uxy,uxz]=gradient(data(:,:,:,1),voxsz(2)/maxvoxsz,voxsz(1)/maxvoxsz,voxsz(3)/maxvoxsz); % X displacement map
% [uyx,uyy,uyz]=gradient(data(:,:,:,2),voxsz(2)/maxvoxsz,voxsz(1)/maxvoxsz,voxsz(3)/maxvoxsz); % Y displacement map
% [uzx,uzy,uzz]=gradient(data(:,:,:,3),voxsz(2)/maxvoxsz,voxsz(1)/maxvoxsz,voxsz(3)/maxvoxsz); % Z displacement map
% found this at
% www.mathworks.com/matlabcentral/newsreader/view_thread/319522
Jacob=(1+uxx).*((1+uyy).*(1+uzz)-(uyz.*uzy))-uxy.*((uyx.*(1+uzz))-(uyz.*uzx))+uxz.*((uyx.*uzy)-((1+uyy).*uzx));

% Jacobian need to be flipped in dimension X (2) and Z (3) to match
% exp_label

% Jacob=flipdim(Jacob,3); % flip Z dimension
% Jacob=flipdim(Jacob,2); % flip X dimension

% another solution based on Reinhardt et al., Med Image Anal. 2008;
% 12(6):752-63 and Wikipedia: Determinant
% stats generated slighty different results between methods (i.e. Jacob1
% and Jacob)
% aei=(1+uxx).*(1+uyy).*(1+uzz);
% bfg=(uxy.*uyz.*uzx);
% cdh=(uxz.*uyx.*uzy);
% ceg=(uxz.*(1+uyy).*uzx);
% bdi=(uxy.*uyx.*(1+uzz));
% afh=((1+uxx).*uyz.*uzy);
% Jacob=(aei+bfg+cdh)-(ceg+bdi+afh);

% Calculate stats of Jacobian
vJ=Jacob(mask);
Jstats=[mean(vJ), median(vJ), std(vJ) skewness(vJ) kurtosis(vJ)];assignin('base','Jstats',Jstats);

% Save Jacobian
cmi_save(0,Jacob,{'Jacob'},voxsz.*matrix);


toc;
