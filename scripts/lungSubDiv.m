function L = lungSubDiv(obj)
% Divides input mask or label matrix into three equally z-spaced sections
% 	for lung segmentation

iflag = 0;
if islogical(obj)
    disp('Finding connected regions in binary mask ...');
    mask = bwlabeln(obj);
elseif isa(obj,'CMIclass') && obj.img.mask.check
    iflag = 2;
    mask = obj.img.mask.mat;
elseif isa(obj,'ImageClass') && obj.mask.check
    iflag = 1;
    mask = obj.mask.mat;
end

ni = max(mask(:));
if ni==1
    disp('Dividing entire mask into 3 z-sections.');
elseif ni==2
    disp('Dividing each lung into 3 z-sections.');
else
    warning('More than 2 regions found ...');
end

d = size(mask);
L = zeros(d);
for i = 1:ni
    lmx = 3*(i-1);
    % Isolate single region
        tmask = mask==i;
    % z-index vector
        iz = squeeze(max(max(tmask,[],1),[],2));
    % calculate separation coordinates
        sep = linspace(find(iz,1,'first'),find(iz,1,'last'),4);
        fprintf('Region %u: %4.1f  %4.1f  %4.1f  %4.1f\n',i,sep);
    % Initialize temp label matrix
        dmask = double(tmask);
    % First section
        iz = (1:d(3))<sep(2);
        dmask(:,:,iz) = tmask(:,:,iz)*(lmx+1);
    % Second section
        iz = ((1:d(3))>=sep(2)) & ((1:d(3))<sep(3));
        dmask(:,:,iz) = tmask(:,:,iz)*(lmx+2);
    % Third section
        iz = (1:d(3))>=sep(3);
        dmask(:,:,iz) = tmask(:,:,iz)*(lmx+3);
    % Add to label matrix:
        L(tmask) = dmask(tmask);
end

if iflag==1 % ImageClass object
    obj.imgAppend(L,{'LungSubDiv'});
elseif iflag==2 % CMIclass object
    obj.imgAppend(L,{'LungSubDiv'});
end