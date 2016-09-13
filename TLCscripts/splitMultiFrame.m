function splitMultiFrame20140819(sourcefile, destidir)
% called by copyMultiFrame()
% split enhanced dicom multiframe image files into single frame image files
% 9/21/06   HG
% 6/23/2012 TLC Bugfix if no series description
% 20140211  TLC Some changes to keep desired fields discovered missing when reading enhanced dicom.
% MD: 20140605 -- safety-checks for remove-fields
% TLC: 20140819 -- Deal with the occasional SecondaryCapture


%disp(' ');
%disp(['Working on destination ' destidir]);
warning off Images:genericDICOM;
    
    series_data = dicomread(sourcefile);        
	series_dcminfo = dicominfo(sourcefile);
    
    %save('C:\tlchenev_stuff\lastser.mat', 'series_dcminfo');
    if (isfield(series_dcminfo, 'SeriesDescription'))
        serdisc = series_dcminfo.SeriesDescription;
    else
        serdisc = 'No Description';
    end
    if (isfield(series_dcminfo, 'SeriesNumber'))
        sernumb = series_dcminfo.SeriesNumber;
        disp(['Series # ' num2str(sernumb) '; Description: ' serdisc]);
    end
    
    % The following UM script merely removes several public tags. Some pub/priv removed below too.
    %frame_dcminfo = get_frame_info2013b(series_dcminfo);
    
    % 20140819:  Check for secondary capture
    
    if ( isfield(series_dcminfo,'DateOfSecondaryCapture') )
        % Deal with ScreenCapture:
        SecondaryCapture = 'y';
        disp(['    FYI: Series # ' num2str(sernumb) ' is a SecondaryCapture.']);
        acquisitionDateTime = [series_dcminfo.AcquisitionDate series_dcminfo.StudyTime];
        destifile = [destidir, acquisitionDateTime, '.', sprintf('%05d', 1), '.dcm'];
        tempserdesc = series_dcminfo.SeriesDescription;
        series_dcminfo.SeriesDescription = [tempserdesc '_SecondaryCapture'];
        dicomwrite(series_data, destifile, series_dcminfo, 'WritePrivate', true, 'CreateMode', 'copy');
        
    else
        
        SecondaryCapture = 'n';
        frame_dcminfo = remove_frame_info(series_dcminfo); % 20140216
        %end %if DateOfSecondaryCapture


        %remove trouble tags, otherwise,they cause dicomwrite failure
        if (isfield(frame_dcminfo, 'ProcedureCodeSequence'))
            frame_dcminfo = rmfield(frame_dcminfo, 'ProcedureCodeSequence'); %45
        end
        if (isfield(frame_dcminfo, 'PerformedProtocolCodeSequence'))
            frame_dcminfo = rmfield(frame_dcminfo, 'PerformedProtocolCodeSequence'); %158
        end
        if (isfield(frame_dcminfo,'Private_2001_105f'))
            frame_dcminfo = rmfield(frame_dcminfo, 'Private_2001_105f'); %189
        end

        if (isfield(frame_dcminfo,'Private_2005_1080')) % MD: 20140605 -- safety-checks
            frame_dcminfo = rmfield(frame_dcminfo, 'Private_2005_1080'); %255
        end
        if (isfield(frame_dcminfo,'Private_2005_1083'))
            frame_dcminfo = rmfield(frame_dcminfo, 'Private_2005_1083'); %255
        end
        if (isfield(frame_dcminfo,'Private_2005_1084'))
            frame_dcminfo = rmfield(frame_dcminfo, 'Private_2005_1084'); %257
        end
        if (isfield(frame_dcminfo,'Private_2005_1085'))
            frame_dcminfo = rmfield(frame_dcminfo, 'Private_2005_1085'); %258
        end
        if (isfield(frame_dcminfo,'Private_2005_109e'))
            frame_dcminfo = rmfield(frame_dcminfo, 'Private_2005_109e'); %260
        end
        if (isfield(frame_dcminfo,'Private_2005_1371'))
            frame_dcminfo = rmfield(frame_dcminfo, 'Private_2005_1371'); %310
        end
        if (isfield(frame_dcminfo,'Private_2005_1402'))
            frame_dcminfo = rmfield(frame_dcminfo, 'Private_2005_1402'); %321
        end
        if (isfield(frame_dcminfo,'Private_2005_140e'))
            frame_dcminfo = rmfield(frame_dcminfo, 'Private_2005_140e'); %323
        end
        if (isfield(frame_dcminfo,'Private_2005_140f'))
            frame_dcminfo = rmfield(frame_dcminfo, 'Private_2005_140f'); %324
        end

        % get only the tags which are in the single frame dicominfo 
        % from SharedFunctionalGroupsSequence


        sharedFields = series_dcminfo.SharedFunctionalGroupsSequence.Item_1;

        frame_dcminfo.RepetitionTime = sharedFields.MRTimingAndRelatedParametersSequence.Item_1.RepetitionTime;
        frame_dcminfo.EchoTrainLength = sharedFields.MRTimingAndRelatedParametersSequence.Item_1.EchoTrainLength;
        frame_dcminfo.FlipAngle = sharedFields.MRTimingAndRelatedParametersSequence.Item_1.FlipAngle;
        frame_dcminfo.NumberOfAverages = sharedFields.MRAveragesSequence.Item_1.NumberOfAverages;
        frame_dcminfo.NumberOfPhaseEncodingSteps = sharedFields.Private_2005_140e.Item_1.NumberOfPhaseEncodingSteps;

        % Start 20140211 ...
        if (isfield(sharedFields,'MRFOVGeometrySequence'))
            roworcolumn = sharedFields.MRFOVGeometrySequence.Item_1.InPlanePhaseEncodingDirection;
            frame_dcminfo.InPlanePhaseEncodingDirection = roworcolumn(1:3);
            % frame_dcminfo = setfield(frame_dcminfo,'InPlanePhaseEncodingDirection',roworcolumn(1:3)); % Prefer to keep "COL" in lieu of "COLUMN"
        end

        if (isfield(sharedFields,'MRFOVGeometrySequence')) % PMS only
            if (isfield(sharedFields.MRFOVGeometrySequence,'Item_1'))
                if (isfield(sharedFields.MRFOVGeometrySequence.Item_1,'PercentPhaseFieldOfView'))
                    frame_dcminfo.PercentPhaseFieldOfView = sharedFields.MRFOVGeometrySequence.Item_1.PercentPhaseFieldOfView;
                end % if 'PercentPhaseFieldOfView'
            end % if 'Item_1'
        end % if 'MRFOVGeometrySequence'


        % ... End 20140211


        % get all the frame specific tags in PerFrameFunctionalGroupsSequence
        frameRecords = series_dcminfo.PerFrameFunctionalGroupsSequence;
        frameRecordsFieldNames = fieldnames(frameRecords);
        numOfFrame = size(frameRecordsFieldNames,1);

        % disp(['FYI: Number of frames this series is = ' num2str(numOfFrame)]);

        for i = 1: numOfFrame
            currentFrame = getfield(frameRecords,frameRecordsFieldNames{i});
            currentFrameFieldNames = fieldnames(currentFrame);
            numOfFrameField = size(currentFrameFieldNames, 1);

            for j = 1:numOfFrameField
                currentField = getfield(currentFrame,currentFrameFieldNames{j});
                %temp = currentField
                if (isa(currentField, 'struct'))
                    if (isa(currentField.Item_1, 'struct'))
                        frameTagsNames = fieldnames(currentField.Item_1);
                        numOfFrameTag = size(frameTagsNames,1);

                        for k = 1:numOfFrameTag
                            currentTag = getfield(currentField.Item_1,frameTagsNames{k});
                            frame_dcminfo = setfield(frame_dcminfo, frameTagsNames{k}, currentTag);               
                        end

                     end

                end
            end

            % get frame specific image data
            frame_data = series_data(:,:,:,i);

            % generate single frame image
            acquisitionDateTime = frame_dcminfo.AcquisitionDateTime;
            destifile = [destidir, acquisitionDateTime, '.', sprintf('%05d', i), '.dcm'];
            dicomwrite(frame_data, destifile, frame_dcminfo, 'WritePrivate', true, 'CreateMode', 'copy');

        end % for i
    
    end % if SecondaryCapture
 
  

    
        
  
    