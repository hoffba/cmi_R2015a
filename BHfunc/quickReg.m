function quickReg(cmiRef,cmiHom)


% Determine file names:
bpath = uigetdir(pwd,'Select timepoint directory:');
[~,bname] = fileparts(bpath);
reffname = fullfile(bpath,[bname,'_Exp.mhd']);
homfname = fullfile(bpath,[bname,'_Ins.mhd']);
voifname = fullfile(bpath,[bname,'_Exp_lungsVOI.mhd']);

% Load images and VOI:
disp('Loading Ref ...')
cmiRef.loadImg(false,reffname);
disp('Loading Hom ...')
cmiHom.loadImg(false,homfname);
disp('Loading RefVOI ...')
cmiRef.loadMask(voifname);

% Dilate VOI:
disp('Dilating Mask ...')
cmiRef.img.mask.morph('dilate',[3,3]);

% Start Coregistration:
regCMIElx(cmiRef,cmiHom);