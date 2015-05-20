% Script to adjust gapped data, inserting blank slices
function ind = adjustGap(cmiObj,targetdz,fillval)

if nargin<3
    fillval = 0;
end

dnew = floor(cmiObj.img.voxsz(3)/targetdz);
d = cmiObj.img.dims;
fov = d(1:3).*cmiObj.img.voxsz;
d(3) = dnew * d(3);
% ind = round(dnew/2):dnew:d(3); 

%% CJG added round b/c I got decimals
answer = questdlg('Is this Leuven Data?');

if strcmp(answer,'Yes')
    ind1 = round(dnew/2):dnew*2:d(3);
    ind2 = ind1+1;
    ind = zeros(1,2*numel(ind1));
    ind(1:2:end) = ind1;
    ind(2:2:end) = ind2;
elseif strcmp(answer,'No')
    ind = round(dnew/2):dnew:d(3);
end

%%
if ~strcmp(answer,'Cancel')
    omask = cmiObj.img.mask.mat;
    
    % Set new image:
    iimg = zeros(d)+fillval;
    iimg(:,:,ind,:) = cmiObj.img.mat;
    cmiObj.setImg(iimg,cmiObj.img.labels,fov);
    
    % Set new mask:
    if ~isempty(omask)
        imask = false(d(1:3));
        imask(:,:,ind) = omask;
        cmiObj.img.mask.merge('replace',imask);
    end
    
    disp(['Original slice numbers are now: ',num2str(ind)])
end