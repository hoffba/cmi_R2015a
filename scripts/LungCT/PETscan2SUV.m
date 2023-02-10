function [SUV,info] = PETscan2SUV(dcmfolder)
SUV = [];

% Find DICOM files in folder
fn = dir(fullfile(dcmfolder,'*.dcm'));

% Loop over slices:
slc_pos = [];

for i = 1:numel(fn)
    % Read DICOM
    info(i) = dicominfo(fullfile(dcmfolder,fn(i).name));
    if ~strcmp(info.Modality,'PT')
        return;
    end
    img = double(dicomread(info));
    
    % Find relevant header info:
    t_acq = str2double(info.AcquisitionTime);
    t_rad = str2double(info.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime);
    half_life = info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideHalfLife/ 60; % [min]
    total_dose = info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose;
    
    % Calculate SUV factor
    delta_time = (t_acq - t_rad) / 100; % [min]
    corrected_dose = total_dose * exp(- delta_time * log(2) / half_life); % [Bq] 
    SUV_factor = info.RescaleSlope * info.PatientWeight * 1000 / corrected_dose; % [g/Bq] = [] * [Kg]* 1000[g/kg] / [Bq] 

    % Create SUV image
    SUV = cat(3,SUV,img * SUV_factor); %[g/ml] = [Bq/ml] * [g/Bq]
    slc_pos = [slc_pos;info.SliceLocation];
end

[~,ix] = sort(slc_pos);
SUV = SUV(:,:,ix);




