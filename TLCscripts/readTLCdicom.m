function [img,fov,label] = readTLCdicom(cmiObj)
% Load DICOM files from the clinical research MRI (Philips)

[imdat,~,~,~,~,fov(1),fov(2),fov(3),~,~,unique_te,~,~,~,~] = ...
                                                readdicom7_combo(0,1,'fp');

img = getsafield(imdat);

% Determine labels (index+1):
icode = {'Unknown', 'MAG',      'RE',       'IM',       'Phase',...
         'ADC',     'FA',       'Vel',      'Water',    'InPhase',...
         'OutPhase','Fat',      'T1',       'T2',       'Other',...
         'ScrnShot','Secondary','FatFract', 'T2star'};
imgtyp = getsafield(imdat,'imgtyp');
imgtyp = imgtyp(1,:);
label = icode(imgtyp+1);

mte = unique_te(unique_te>0);
disp(['TE =',sprintf('  %.3f',mte)]);
assignin('base','te',mte);

if nargin && isa(cmiObj,'CMIclass')
    cmiObj.setImg(img,label,double(fov));
end