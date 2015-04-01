function [ Yc_vol , Yc , Ym ] = vol_homocor_cjg_3D(img) 
% -----------------------------------------------------------
%	to correct for image nonuniformity
%	first fills in holes, then generates a low frequency 
%	background image and makes correction factor
%	input image is stored in img, output is Iout
%
%	homocor.m  rev 0    4/8/00	          Gary Glover
%       ----------------------------------------------------- 
%
%   Modified to deal with SPM99 analize format image files.
%   Produces a 'homocor.img' volume in the working directory.
%
%   Works with volumes acquired in any orientation, The inplane
%   is taken to be the plane where two of the volume dimensions
%   match. If all dimensions match, it will default to axial.
% -----------------------------------------------------------
%   @(#)vol_homocor.m    1.10  Kalina Christoff    2000-05-29
% -----------------------------------------------------------
%   Notes from Craig Galban 17June2005 ----------------------
%   From trial5 scan 2 proc 3, <10% difference from original to corrected
%   was obtained when using thres=5, fspecial('gaussian',21,Ist) where Ist
%   is std(std(std(img)))*10.

% get the path
% ------------
% vol_path=spm_get(1, 'img', 'Please select image to correct');
% [vol_file,vol_path]=uigetfile('*.*','Please select image to correct');
% load the volume
% ---------------
% if ischar(vol_path)
%     vol=load([vol_path,vol_file]);
%     %   vol = spm_vol(vol_path);
% end
% make a guess as to what orientation are the slices
% eg. vol.dim=[256 256 124 4] probably axial, etc
% --------------------------------------------------
vol.dim=size(vol.eye_data);
if       vol.dim(1)>=1 % probably axial
            inpl_crd=[1 2 3]; thrpl_crd=3;
            orientation='3D';

%  elseif  vol.dim(2)==vol.dim(3) % probably sagittal
%             inpl_crd=[2 3]; thrpl_crd=1;
%             orientation='sagittal';
% 
%  elseif  vol.dim(1)==vol.dim(3) % probably coronal
%             inpl_crd=[1 3]; thrpl_crd=2;
%             orientation='coronal'; 
end


  fprintf('\nAssuming the planes were acquired in ');
  fprintf(orientation); fprintf( ' orientation.\n\n');
      
foo = vol.dim(inpl_crd);
npix = foo(1);
np2 = npix/2; % one half dimension

%thres = input('gimme threshold % [5] =');
%if(isempty(thres))
% thres = 5;
%end;
thres = handles.thres*.01;

% Y = spm_read_vols(vol,0);
Y = vol.eye_data;

% initialize the corrected output volume Yc
Yc=zeros(vol.dim(1:3));

fprintf('Please wait - now working on plane    ');

% begin looping through slices
% ----------------------------
% for slice = 1:vol.dim(thrpl_crd);
     
%     if slice<10; fprintf('\b%d',slice); 
%      elseif slice<100,  fprintf('\b\b%d',slice); 
%      elseif slice<1000, fprintf('\b\b\b%d',slice);
%     end

    % surprisingly, Y needs to be squeezed to remove the singleton 
    % dimensions for coronal or sagittal... interesting command.
    % ---------------------------------------------------------------
    if thrpl_crd==1; img=squeeze(Y(slice,:,:)); end;  
    if thrpl_crd==2; img=squeeze(Y(:,slice,:)); end;
    if thrpl_crd==3; img=squeeze(Y(:,:,:)); end; % though no need here

    Imax = max(max(max(img))); % calculate max
    Ithr = Imax*thres; % multiply max with threshold
    
    %  fill in holes

    mask = (img>Ithr); % create mask
    Imask = img.*mask;
    Iave = sum(sum(sum(Imask)))/sum(sum(sum(mask))); % weighted average
    Ifill = Imask + (1-mask).*Iave; % fill all points that are zero in the mask add Iave
    Ist=std(std(std(img)))*10;
    %  make low freq image

%     z = fftshift(fftn(Ifill));  % take fft2

    %ndef = 3;
    ndef = np2/8;
    n = ndef;

    y2 = (1:np2).^2;  % use 15 as default win for guassian
    a2 = 1/n;
    
    fil = fspecial('gaussian',21,Ist);
    Ilow=imfilter(Ifill,fil);

%     Ilow = abs(ifftn(fftshift(zl)));

    %corrected image

    Ilave = sum(sum(sum(Ilow.*mask)))/sum(sum(sum(mask)));
    Icorf = (Ilave./Ilow).*mask;  % correction factor
    Iout = img.*Icorf;
    Ym=Y.*mask;

    % update the volume with the corrected values
    % -------------------------------------------
    if thrpl_crd==1; Yc(slice,:,:) = Iout; end
    if thrpl_crd==2; Yc(:,slice,:) = Iout; end
    if thrpl_crd==3; Yc(:,:,:) = Iout; end
  
    clear img;
  
% end looping through slices
% ---------------------------
% end

    fprintf('\n\nDone.\n\n')

% write corrected image (homocor.img)
% ==================================

Yc_vol = vol;
Yc_vol.fname   ='homocor';
Yc_vol.descrip ='homocor';

% spm_write_vol(Yc_vol,Yc);