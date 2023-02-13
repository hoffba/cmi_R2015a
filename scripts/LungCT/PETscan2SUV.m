function [SUV,fov,orient,info_out] = PETscan2SUV(dcmfolder)
SUV = [];
info_out = [];

% Find DICOM files in folder
fn = dir(fullfile(dcmfolder,'*.dcm'));

% Loop over slices:
ns = numel(fn);
slc_pos = nan(ns,1);
for i = 1:ns
    % Read DICOM
    info(i) = dicominfo(fullfile(dcmfolder,fn(i).name));
    if ~strcmp(info(i).Modality,'PT')
        return;
    end
    img = double(dicomread(info(i)));
    
    % Find relevant header info:
    t_acq = str2double(info(i).AcquisitionTime);
    t_rad = str2double(info(i).RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime);
    half_life = info(i).RadiopharmaceuticalInformationSequence.Item_1.RadionuclideHalfLife/ 60; % [min]
    total_dose = info(i).RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose;
    
    % Calculate SUV factor
    delta_time = (t_acq - t_rad) / 100; % [min]
    corrected_dose = total_dose * exp(- delta_time * log(2) / half_life); % [Bq] 
    SUV_factor = info(i).RescaleSlope * info(i).PatientWeight * 1000 / corrected_dose; % [g/Bq] = [] * [Kg]* 1000[g/kg] / [Bq] 

    % Create SUV image
    SUV = cat(3,SUV,img * SUV_factor); %[g/ml] = [Bq/ml] * [g/Bq]
    slc_pos(i) = info(i).SliceLocation;
end

[slc_pos,ix] = sort(slc_pos);
info = info(ix);
info_out = info(1);
SUV = SUV(:,:,ix);
fov = (slc_pos(end) - slc_pos(1)) * (1 + 1/(ns-1));
orient = [ [ info_out.ImageOrientationPatient(1:3)*info_out.PixelSpacing(1) ,...
             info_out.ImageOrientationPatient(4:6)*info_out.PixelSpacing(2) ,...
             (info(end).ImagePositionPatient - info(1).ImagePositionPatient)/(ns-1) ,...
             info_out.ImagePositionPatient ] ;...
             0 0 0 1];





