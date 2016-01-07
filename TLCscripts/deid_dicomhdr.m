function di1_lite = deid_dicomhdr(di1);  
% Extract some key technical factors from dicom header, toss the rest.
% Done to remove PHI
% Thomas L Chenevert UMich Aug 2011.
% TLC July 2012.  Add series time
% TLC 20141125.  ReconstructionDiameter gone in R5.  Use Rows&Columns*PixelSpacing


% Create an abreviated dicom header.  Keep non-PHI and only several key parameters.
di1_lite.ImageType = getexistfield(di1, 'ImageType');
di1_lite.Manufacturer = getexistfield(di1, 'Manufacturer');
di1_lite.StationName = getexistfield(di1, 'StationName');
di1_lite.StudyDescription = getexistfield(di1, 'StudyDescription');
di1_lite.ScanningSequence = getexistfield(di1, 'ScanningSequence');
di1_lite.MRAcquisitionType = getexistfield(di1, 'MRAcquisitionType');
di1_lite.SliceThickness = getexistfield(di1, 'SliceThickness');
di1_lite.RepetitionTime = getexistfield(di1, 'RepetitionTime');
di1_lite.EchoTime = getexistfield(di1, 'EchoTime');
di1_lite.NumberOfAverages = getexistfield(di1, 'NumberOfAverages');
di1_lite.ImagingFrequency = getexistfield(di1, 'ImagingFrequency');
di1_lite.SpacingBetweenSlices = getexistfield(di1, 'SpacingBetweenSlices');
di1_lite.NumberOfPhaseEncodingSteps = getexistfield(di1, 'NumberOfPhaseEncodingSteps');
di1_lite.PercentSampling = getexistfield(di1, 'PercentSampling');
di1_lite.PercentPhaseFieldOfView = getexistfield(di1, 'PercentPhaseFieldOfView');
di1_lite.PixelBandwidth = getexistfield(di1, 'PixelBandwidth');
di1_lite.DeviceSerialNumber = getexistfield(di1, 'DeviceSerialNumber');
di1_lite.SoftwareVersion = getexistfield(di1, 'SoftwareVersion');
di1_lite.ReconstructionDiameter = getexistfield(di1, 'ReconstructionDiameter');
di1_lite.AcquisitionMatrix = getexistfield(di1, 'AcquisitionMatrix');
di1_lite.FlipAngle = getexistfield(di1, 'FlipAngle');
di1_lite.NumberOfTemporalPositions = getexistfield(di1, 'NumberOfTemporalPositions');
di1_lite.Rows = getexistfield(di1, 'Rows');
di1_lite.Columns = getexistfield(di1, 'Columns');
di1_lite.PixelSpacing = getexistfield(di1, 'PixelSpacing');
di1_lite.PerformedStationAETitle = getexistfield(di1, 'PerformedStationAETitle');
di1_lite.SeriesTime = getexistfield(di1, 'SeriesTime');
ps = di1_lite.PixelSpacing;
di1_lite.fovrows = (di1_lite.Rows) * ps(1);
di1_lite.fovcols = (di1_lite.Rows) * ps(2);
