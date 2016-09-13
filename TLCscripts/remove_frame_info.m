function frame_info = remove_frame_info(series_info)
% Used to be called "get_frame_info" and was called by splitMultiFrame()
% but since purpose of this script is to remove fields anticipated to be
% problematic as we convert enhanced dicom to non-enhanced dicom, we
% decided to rename it to "remove_frame_info" and only keep calls to
% "rmfield" that either cause problems, or are not currently needed.
% DO NOT delete the line, but comment-out but leave visible for future reference.
% Circa install of 2013b (20140211) we could see many useful fields that
% were being removed previously.  Put those field (and others) back in
% since down-stream script can use those tags.
% Original intent of HG: remove dicom tags which are in the enhanced dicominfo('multi-frame image file'), 
% but not in the standard dicominfo('single frame image file')
% 9/21/06   HG
% 10/12/2007 TLC Check if certain fields exist before removing them, else
% an error.  See, for example, changes below for field 'ImageComments'.
% TLC Mar 04 , 2011 check if 'VolumetricProperties' isfield
% TLC Aug12, 2011 check if 'RectilinearPhaseEncodeReordering' isfield
% TLC 20140211. Circa 2013b install. Blot-out to NOT remove some fields.
% Note, the more fields retained (i.e. rmfield is commented-out), the slower the script will run.

    frame_info = series_info;
    frame_info = rmfield(frame_info, 'CreatorVersionUID');
    frame_info = rmfield(frame_info, 'PixelPresentation');
    if (isfield(frame_info, 'VolumetricProperties'))
        frame_info = rmfield(frame_info, 'VolumetricProperties');
    end
    
    if (isfield(frame_info, 'VolumeBasedCalculationTechnique'))
        frame_info = rmfield(frame_info, 'VolumeBasedCalculationTechnique');
    end
    
    frame_info = rmfield(frame_info, 'ComplexImageComponent');
    frame_info = rmfield(frame_info, 'AcquisitionContrast');
    frame_info = rmfield(frame_info, 'ContentQualification');
    frame_info = rmfield(frame_info, 'PulseSequenceName');
    frame_info = rmfield(frame_info, 'EchoPulseSequence');
    frame_info = rmfield(frame_info, 'MultiplanarExcitation');
    frame_info = rmfield(frame_info, 'PhaseContrast');
    frame_info = rmfield(frame_info, 'TimeOfFlightContrast');
    frame_info = rmfield(frame_info, 'Spoiling');
    frame_info = rmfield(frame_info, 'SteadyStatePulseSequence');
    frame_info = rmfield(frame_info, 'EchoPlanarPulseSequence');
    frame_info = rmfield(frame_info, 'MagnetizationTransfer');
    frame_info = rmfield(frame_info, 'T2Preparation');
    frame_info = rmfield(frame_info, 'BloodSignalNulling');
    frame_info = rmfield(frame_info, 'SaturationRecovery');
    frame_info = rmfield(frame_info, 'SpectrallySelectedSuppression');
    frame_info = rmfield(frame_info, 'SpectrallySelectedExcitation');
    frame_info = rmfield(frame_info, 'SpatialPresaturation');
    frame_info = rmfield(frame_info, 'Tagging');
    frame_info = rmfield(frame_info, 'OversamplingPhase');
    frame_info = rmfield(frame_info, 'GeometryOfKSpaceTraversal');
    frame_info = rmfield(frame_info, 'SegmentedKSpaceTraversal');
     
%     %frame_info = rmfield(frame_info, 'RectilinearPhaseEncodeReordering');
    if (isfield(frame_info, 'RectilinearPhaseEncodeReordering'))
        frame_info = rmfield(frame_info, 'RectilinearPhaseEncodeReordering');
    end
  
    frame_info = rmfield(frame_info, 'TagThickness');
%     if (isfield(frame_info, 'PartialFourierDirection')) % Keep 20140211
%         frame_info = rmfield(frame_info, 'PartialFourierDirection');
%     end
    frame_info = rmfield(frame_info, 'CardiacSynchronizationTechnique');
    frame_info = rmfield(frame_info, 'TransmitCoilType');
    frame_info = rmfield(frame_info, 'ChemicalShiftReference');
    frame_info = rmfield(frame_info, 'MRAcquisitionFrequencyEncodingSteps');
    frame_info = rmfield(frame_info, 'Decoupling');
    frame_info = rmfield(frame_info, 'KSpaceFiltering');
%     % frame_info = rmfield(frame_info, 'ParallelReductionFactorInPlane'); % Keep 20140211.
%     frame_info = rmfield(frame_info, 'AcquisitionDuration');      % Keep 20140211
%     % frame_info = rmfield(frame_info, 'ParallelAcquisition');  % Keep 20140211.
%     if (isfield(frame_info, 'ParallelAcquisitionTechnique')) % Keep 20140211.
%         frame_info = rmfield(frame_info, 'ParallelAcquisitionTechnique');
%     end
%     frame_info = rmfield(frame_info, 'PartialFourier'); % Keep 20140211
    frame_info = rmfield(frame_info, 'VelocityEncodingDirection');
    frame_info = rmfield(frame_info, 'VelocityEncodingMinimumValue');
    frame_info = rmfield(frame_info, 'NumberOfKSpaceTrajectories');
    if (isfield(frame_info, 'CoverageOfKSpace'))
        frame_info = rmfield(frame_info, 'CoverageOfKSpace');
    end
    frame_info = rmfield(frame_info, 'ResonantNucleus');
    frame_info = rmfield(frame_info, 'FrequencyCorrection');
    frame_info = rmfield(frame_info, 'ParallelReductionFactorOutOfPlane');
    frame_info = rmfield(frame_info, 'ParallelReductionFactorSecondInPlane');
    frame_info = rmfield(frame_info, 'RespiratoryMotionCompensationTechnique');
    frame_info = rmfield(frame_info, 'BulkMotionCompensationTechnique');
    frame_info = rmfield(frame_info, 'ApplicableSafetyStandardAgency');
%     frame_info = rmfield(frame_info, 'SpecificAbsorptionRateDefinition'); % Keep 20140211.
    frame_info = rmfield(frame_info, 'GradientOutputType');
%     frame_info = rmfield(frame_info, 'SpecificAbsorptionRateValue'); % Keep 20140211.
    frame_info = rmfield(frame_info, 'GradientOutput');
    frame_info = rmfield(frame_info, 'WaterReferencedPhaseCorrection');
    frame_info = rmfield(frame_info, 'MRSpectroscopyAcquisitionType');
%     frame_info = rmfield(frame_info, 'MRAcquisitionPhaseEncodingStepsInPlane'); % Keep 20140211.
    frame_info = rmfield(frame_info, 'RFEchoTrainLength');
    frame_info = rmfield(frame_info, 'GradientEchoTrainLength');
    
    %frame_info = rmfield(frame_info, 'ImageComments');
    if (isfield(frame_info, 'ImageComments'))
        frame_info = rmfield(frame_info, 'ImageComments');
    end
    frame_info = rmfield(frame_info, 'FrameLaterality');
    frame_info = rmfield(frame_info, 'DimensionOrganizationSequence');
    frame_info = rmfield(frame_info, 'DimensionIndexSequence');
    frame_info = rmfield(frame_info, 'NumberOfFrames');
    frame_info = rmfield(frame_info, 'BurnedInAnnotation');
    if (isfield(frame_info, 'LUTExplanation'))
        frame_info = rmfield(frame_info, 'LUTExplanation');
    end
    frame_info = rmfield(frame_info, 'DataPointRows');
    frame_info = rmfield(frame_info, 'DataPointColumns');
    frame_info = rmfield(frame_info, 'AcquisitionContextSequence');
    frame_info = rmfield(frame_info, 'SharedFunctionalGroupsSequence');
    frame_info = rmfield(frame_info, 'PerFrameFunctionalGroupsSequence');
    