% ImageIO method
function cmiInfo = parseInfo(self,info,itype)
% Input raw image info from loader functions
% Returns structure containing CMI-required information

cmiInfo = self.definfo;
switch itype
    case 1 % Analyze
    case 2 % MetaIO
        if isfield(info,'x_Institution')
            cmiInfo.Institution = info.x_Institution;
        end
        if isfield(info,'x_Manufacturer')
            cmiInfo.Manufacturer = info.x_Manufacturer;
        end
        if isfield(info,'x_Model')
            cmiInfo.Model = info.x_Model;
        end
        if isfield(info,'x_FieldStrength')
            cmiInfo.FieldStrength = info.x_FieldStrength;
        end
        if isfield(info,'x_StationName')
            cmiInfo.StationName = info.x_StationName;
        end
        if isfield(info,'x_SerialNo')
            cmiInfo.SerialNo = info.x_SerialNo;
        end
        if isfield(info,'x_SoftwareVersion')
            cmiInfo.SoftwareVer = info.x_SoftwareVersion;
        end
        if isfield(info,'x_StudyDescription')
            cmiInfo.StudyDescription = info.x_StudyDescription;
        end
        if isfield(info,'x_SeriesDescription')
            cmiInfo.SeriesDescription = info.x_SeriesDescription;
        end
        if isfield(info,'x_ImageTypeString')
            cmiInfo.ImageType = info.x_ImageTypeString;
        end
        if isfield(info,'x_ExternalComment_opt_')
            cmiInfo.Comment = info.x_ExternalComment_opt_;
        end
        if isfield(info,'StationName')
            cmiInfo.Modality = info.StationName;
        end
        cmiInfo.Dim(1:length(info.DimSize)) = info.DimSize;
        if isfield(info,'x_4thDimLabel') % Tom's custom tags
            lstr = {info.x_4thDimLabel};
        else
            [~,lstr,~] = fileparts(info.FileName);
            lstr = {lstr};
        end
        if length(lstr)~=cmiInfo.Dim(4)
            lstr = strcat(lstr{1},'_',...
                cellfun(@num2str,num2cell(1:cmiInfo.Dim(4)),...
                'UniformOutput',false));
        end
        cmiInfo.Label = lstr;
        if isfield(info,'ElementSpacing')
            cmiInfo.VoxelSpacing(1:length(info.ElementSpacing)) = info.ElementSpacing;
        elseif isfield(info,'ElementSize')
            cmiInfo.VoxelSpacing(1:length(info.ElementSize)) = info.ElementSize;
        end
        if isfield(info,'ElementSize')
            cmiInfo.VoxelSize(1:length(info.ElementSize)) = info.ElementSize;
        end
        if isfield(info,'Offset')
            cmiInfo.Origin = info.Offset;
        end
        if isfield(info,'TransformMatrix')
            cmiInfo.Orientation = reshape(info.TransformMatrix,[3,3]);
        end
        if isfield(info,'AnatomicalOrientation')
            cmiInfo.AnatomicalOrientation = info.AnatomicalOrientation;
        end
    case 3 % AVS
    case 4 % FDF
    case 5 % VFF
    case 6 % DICOM
        if isfield(info,'ClinicalTrialSiteName')
            cmiInfo.Institution = info.ClinicalTrialSiteName;
        end
        if isfield(info,'Manufacturer')
            cmiInfo.Manufacturer = info.Manufacturer;
        end
        if isfield(info,'ManufacturerModelName')
            cmiInfo.Model = info.ManufacturerModelName;
        end
        if isfield(info,'MagneticFieldStrength')
            cmiInfo.FieldStrength = info.MagneticFieldStrength;
        end
        if isfield(info,'StationName')
            cmiInfo.StationName = info.StationName;
        end
        if isfield(info,'DeviceSerialNumber')
            cmiInfo.SerialNo = info.DeviceSerialNumber;
        end
        if isfield(info,'SoftwareVersion')
            cmiInfo.SoftwareVer = info.SoftwareVersion;
        end
        if isfield(info,'RequestedProcedureDescription')
            cmiInfo.Comment = info.RequestedProcedureDescription;
        end
        if isfield(info,'Modality')
        	cmiInfo.Modality
        end
        if isfield(info,'StudyDescription')
            cmiInfo.StudyDescription = info.StudyDescription;
        end
        if isfield(info,'SeriesDescription')
            cmiInfo.SeriesDescription = info.SeriesDescription;
        end
        if isfield(info,'ImageType')
            cmiInfo.ImageType = info.ImageType;
        elseif isfield(info,'ImageTypeText')
            cmiInfo.ImageType = info.ImageTypeText;
        end
        if isfield(info,'Rows')
            cmiInfo.Dim(1) = info.Rows;
        end
        if isfield(info,'Columns')
            cmiInfo.Dim(2) = info.Columns;
        end
        if isfield(info,'NumberOfSlices')
            cmiInfo.Dim(3) = info.NumberOfSlices;
        end
        if isfield(info,'ArrayDim')
            cmiInfo.Dim(4) = info.ArrayDim;
        end
        if isfield(info,'PixelSpacing')
            cmiInfo.VoxelSpacing([2,1]) = info.PixelSpacing;
            cmiInfo.VoxelSize([2,1]) = info.PixelSpacing;
        end
        if isfield(info,'SliceThickness')
            cmiInfo.VoxelSize(3) = info.SliceThickness;
        end
        if isfield(info,'SpacingBetweenSlices')
            cmiInfo.VoxelSpacing(3) = info.SpacingBetweenSlices;
        end
        if isfield(info,'ImageOrientationPatient')
            cmiInfo.Orientation = [info.ImageOrientationPatient(1:3),...
                                   info.ImageOrientationPatient(4:6),...
                                   cross(info.ImageOrientationPatient(1:3),...
                                         info.ImageOrientationPatient(4:6))];
        end
        cmiInfo.Origin = info.ImagePositionPatient;
        cmiInfo.Label = info.Label;
        % Need to figure out what DICOM tag correspond to this:
%         if isfield(info,'PatientPosition')
%             cmiInfo.AnatomicalOrientation = info.PatientPosition;
%         end

    case 7 % NifTi
    case 8 % Bruker
    case 9 % MRSolutions
    case 10 % Mask
    case 11 % TIFF
    case 12 % JPEG
    case 13 % Matlab
    case 14 % FID
end

end
