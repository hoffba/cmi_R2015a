function Iout = vol_homocor_bah_3D(Iin,thres) 

if nargin
    % Determine threshold
    if (nargin==1)
        thres = 0.02; % [%]
    elseif thres>1
        disp('Threshold too high!')
    elseif thres<0
        disp('Threshold can''t be negative!')
    end

    % Create thresholded/filtered image
    mask = (Iin > (max(Iin(:))*thres)); % create threshold mask
    Imask = Iin.*mask;
    % fill all points outside the mask with image average:
    Ifill = Imask + (1-mask).*sum(Imask(:))/nnz(mask);
    Ilow = imfilter(Ifill,fspecial('gaussian',21,std(Iin(:))*10));

    % Correct input image:
    tavg = sum(Ilow(mask))/nnz(mask);
    Iout = Iin./Ilow.*mask*tavg;
end
