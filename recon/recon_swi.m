function [MagImg PhaseImg SWImg] = recon_swi(fdir,filt)
% Reconstructs raw SWI MRI images
% Inputs:   fname = location of raw data (includes directories)
%           filt = optional inputs
% Outputs:  MagImg = typical magnitude image
%           HPPhaseImg = high-pass phase image
%           SWImg = susceptibility-weighted image

status = false;
if (nargin == 0)
    % Select raw MRI data
    fdir = uigetdir(pwd,'Select FID to load:');
    if fdir
        status = true;
    end
elseif strcmp(fdir((end-3):end),'.fid')
    status = true;
end


if status
    if (nargin < 2)
        filt = 3;
    end
    nfilt = 4; % power of phase mask
    [kimg,fov,~] = cmi_readFID(fdir);
    if (size(kimg,5) == 1) % ne == 1
        [d(1) d(2) d(3) d(4)] = size(kimg);
        % Determine k-space filtering for noise reduction
        kfilt = zeros(d(1),d(2)); % 2D k-space filter
        switch filt
            case 0
                filt = ones(1,d(1));
            case 1
                filt = ramp(d(1),45,1,30,0,0);
            case 2
                filt = ramp(d(1),35,1,30,0,0);
            otherwise
                filt = ramp(d(1),20,1,10,0,0);
        end
        kx0 = round(d(1)/2); % find center of image
        ky0 = round(d(2)/2);
        for ikx=1:d(1)
            for iky = 1:d(2)
                % Add 1 to avoid index rho of zero.
                rho = round(sqrt((ikx-kx0)^2 + (iky-ky0)^2) + 1);
                kfilt(ikx,iky) = filt(rho);
            end % for iky
        end % for ikx

        xx = zeros(d(1),d(2),d(3),d(4));
        for iarray = 1:d(4)
            for islc = 1:d(3)
                timg = kimg(:,:,islc,iarray);
                ssp0 = fftshift(fft2(timg));
                cref0 = fftshift(fft2(timg.*kfilt));
                xx(:,:,islc,iarray) = (ssp0 .* conj(cref0)) ./ (abs(cref0) + eps);
            end % slci; loop over slices
        end % iarray
        MagImg = abs(xx);
        PhaseImg = angle(xx);
        
        % Calculate phase mask and resulting SWI
        phmask = (PhaseImg/pi + 1);
        phmask(phmask>1) = 1;
        SWImg = MagImg .* (phmask .^ nfilt);

        % Save images as FLD
        outimg = cat(4,MagImg,PhaseImg*10^4,SWImg);
        label = {'Magnitude','Phase(e4)',['SWIp' num2str(nfilt)]};
        cmi_save(false,outimg,label,fov);
    end
end