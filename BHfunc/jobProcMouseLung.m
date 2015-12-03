function jobProcMouseLung(img,fov,label,fname,svdir)
% Processes mouse lung CT images from the Bruker SkyScan 1176
% Performs the following steps:
%   - NOVA filter (de-noise)
%   - Save filtered images as MHD

if ~exist(svdir,'dir')
    disp(svdir)
    error('Save directory does not exist! ');
elseif size(img,4)>1
    error('jobProcMouseLung : img must be 3D! ');
else
    
    % Perform a NOVA filter
    n = 3;
    run = 2;
    sDev0 = 200;
    p = 0;
    d = 3;
    img = bFiltNOVA(img,n,run,sDev0,p,d);
    
    % Save filtered image as MHD:
    saveMHD(fullfile(svdir,fname),img,label,fov);
    
end
