function dcmPatient = read_enh(diropt,selectPatient)
% TLChenever Mar 2011.  Read DVD's list of patients for user to select
% which patient to restore from DVD onto local disk.  Delete MatLab
% image_order.mat residues.
% TLC % Aug 11, 2011:  % Hardcode responses for read_tlcs_dvd(3).  Search "Hardcode Aug 11, 2011".
% TLC 20130930: Create 2nd optional input argument to select which patient to extract.
% TLC 20140211: Try new splitMultiFrame2013b to deal with Enh Dicom issues.
% TLC 20140216: Previously called "read_tlcs_dvd2013b".  Cleaned-up,
% appropriately renamed to "read_enh".
% TLC 20140218: Re-check proper image scaling using read_enh.  Discovered that
% single-image series present img scaling errors, and other gaps in dicom
% header. Improvement, though by forcing all thru "splitMultiFrame" even
% when imgNum = 1.  Therefore, blot-out prior code related to "if imgNum=1".
% TLC 20140808.  Add new "diropt=4" option so this script can be called
%                from "parse_series.m" and the desitnation directory path is shorter.
%                Also, disable all switches based on license. Use 'TOMPC4' settings henceforth.


global lastdcdir;

% ***************************************************************************.
% ***** In case of multi-patient DICOMDIR, Hardcode DesiredPatient Here *****.
%selectPatient = 1; % Hardcode Aug 11, 2011. Or comment-out and use prompt below.
%selectPatient = input('Provide index of patient you want to read from enh DICOM (eg, 1, 2, ... ) ');
% ***************************************************************************.


if (nargin == 0)       % No input argument,
    diropt = 1;        % The default.  May want to add options in the future.
    selectPatient = 1; % Assume you want the 1st patient.
elseif (nargin == 1)   % then diropt is defined by 1st argument
    selectPatient = 1; % Assume you want the 1st patient.
elseif (nargin == 2)
    disp('Two arguments being used as diropt and selectPatient ...');
else
    disp('Too many arguments! Will only use 1st two as diropt and selectPatient');
end % if nargin

disp(['FYI: Using diropt = ' num2str(diropt)]);
disp(['FYI: Using selectPatient = ' num2str(selectPatient)]);


doscaninfo = 'y'; % Hardcode Aug 11, 2011
%doscaninfo = input('Create a scan_info.txt file (y/n) ? ','s'); % Hardcode Aug 11, 2011.

if(doscaninfo == 'y')
    keepmats = 'y'; % Hardcode Aug 11, 2011
    %keepmats = input('Keep image_order.mats (y/n) ? ','s'); % Hardcode Aug 11, 2011
else
    keepmats = 'n';
end % if doscaninfo

tstart = clock; % How long is script?

filecount = 0;
% whichcpu = input('Which CPU are you currently on (1=Fellows_PC; 2=Tom_PC2; 3=TLC_LapTop ) ? ');
% 
% switch whichcpu
%     case 1
%         sourcepath = 'D:\'; % On Research Fellows PC
%         destipath = 'E:\InBox_from_DVD\SimpleDICOM\DicomImages_by_Vendor'; % On Research Fellows PC
%     case 2
%         sourcepath = 'D:\'; % On Tom_PC2
%         %sourcepath = 'D:\Some_Dicom'; % On Tom_PC2
%         %destipath = 'D:\InBox_from_DVD\SimpleDICOM\DicomImages_by_Vendor'; % On Tom_PC2
%         destipath = 'F:\InBox_from_DVD\SimpleDICOM\DicomImages_by_Vendor';% On Tom_PC2
%     case 3
%         sourcepath = 'D:\'; % On TLC_LapTop
%         destipath = 'C:\InBox_from_DVD\SimpleDICOM\DicomImages_by_Vendor'; % On Tom_PC2
% end


% Use these as cheap way to switch-in unconventional, one-off weird character workarounds.
detour1 = 'n';
detour2 = 'n'; % Often 'y' for Siemens data
selectmodality = 'MR';

% %********************************************
% %***** Choose proper "pathology" ************
% %********************************************
% lic = cell(5,1);
% lic(1) = cellstr('498205'); % TOMPC2
% lic(2) = cellstr('498208'); % TOMPC3
% lic(3) = cellstr('498210'); % TOMPC4
% lic(4) = cellstr('498215'); % TLC Lenovo
% lic(5) = cellstr('498216'); % Research Fellows
% lic(6) = cellstr('672490'); % Frank's 64bit Win7 OS system



%thiscpu = char(hostid);
thiscpu = license;

% switch thiscpu
%     
%     case char(lic(1))  % TOMPC2
%         % disp('TOMPC2');
%         if (diropt == 1)
%             % Image Source:
%             sourcepath = 'D:\'; % DVD on Tom_PC2
%             destipath = 'F:\InBox_from_DVD\SimpleDICOM\DicomImages_by_Vendor';% On Tom_PC2
%         else (diropt == 2)
%             disp(['Option ' num2str(diropt) ' is not allowed on this system.  Aborting ... '])
%             return
%         end % if diropt
%         
%     case char(lic(2))  % TOMPC3
%         % disp('TOMPC3');
%         if (diropt == 1)
%             % Image Source: 
%             sourcepath = 'D:\'; % DVD on Tom_PC2
%             destipath = 'E:\InBox_from_DVD\SimpleDICOM\DicomImages_by_Vendor';% On Tom_PC3
%         else (diropt == 2)
%             disp(['Option ' num2str(diropt) ' is not allowed on this system.  Aborting ... '])
%             return
%         end % if diropt
%         
%     case char(lic(3))  % TOMPC4
%         disp('TOMPC4');
        if (diropt == 1)
            % Image Source:
            sourcepath = 'D:\'; % DVD on Tom_PC4
            destipath = 'C:\InBox_from_DVD\SimpleDICOM\DicomImages_by_Vendor';% On Tom_PC4
        elseif (diropt == 2)
            % Manually define Image Source path to place where DICOMDIR lives: 
            sourcepath = 'K:\GABA_FBIRN_FEB_23_2011\MEADOWS_FEB_23_2011'; % DVD on Tom_PC4
            destipath = 'C:\InBox_from_DVD\SimpleDICOM\DicomImages_by_Vendor';% On Tom_PC4
        elseif (diropt == 3)
            % Manually define Image Source path to place where DICOMDIR lives:
            disp('Run This Script In Same Directory as DICOMDIR File. ');
            sourcepath = pwd; % DVD on Tom_PC4
            destipath = pwd;% On Tom_PC4
        elseif (diropt == 4)
            sourcepath = pwd;
            destipath = pwd;
            [basepath,basefolder,baseext] = fileparts(pwd);
            % Define destination exam folder name
            examfolder = [basefolder '_E000'];
        else 
            disp(['Option ' num2str(diropt) ' is not allowed on this system.  Aborting ... ']')
            return
        end % if diropt
%         
%     case char(lic(4))  % TLC Lenovo
%         % disp('Lenovo');
%         if (diropt == 1)
%             % Image Source:
%             sourcepath = 'D:\'; % DVD on Lenovo
%             destipath = 'C:\InBox_from_DVD\SimpleDICOM\DicomImages_by_Vendor';% On Lenovo
%         else
%             disp(['Option ' num2str(diropt) ' is not allowed on this system.  Aborting ... '])
%             return;
%         end % if diropt
%         
%     case char(lic(5))  % Research Fellows
%         % disp('Research Fellows');
%         if (diropt == 1)
%             % Image Source:
%             sourcepath = 'D:\'; % DVD on Tom_PC2
%             destipath = 'C:\InBox_from_DVD\SimpleDICOM\DicomImages_by_Vendor';% On Research Fellows
%         else (diropt == 2)
%             disp(['Option ' num2str(diropt) ' is not allowed on this system.  Aborting ... ']')
%             return
%         end % if diropt
%         
%     case char(lic(6))  % Franks PC
%         disp('Franks PC');
%         if (diropt == 1)
%             % Image Source:
%             sourcepath = 'D:\'; % DVD on Franks PC
%             destipath = 'C:\InBox_from_DVD\SimpleDICOM\DicomImages_by_Vendor';% On Research Fellows
%         elseif (diropt == 2)
%             % Manually define Image Source path to place where DICOMDIR lives:
%             disp('Run This Script In Same Directory as DICOMDIR File. ');
%             sourcepath = pwd; % DVD on Tom_PC4
%             destipath = pwd;  % 
%         else (diropt == 3)
%             disp(['Option ' num2str(diropt) ' is not allowed on this system.  Aborting ... ']')
%             return
%         end % if diropt
%         
% end % switch thiscpu

% disp(['FYI: Source Path is ' sourcepath]);
disp(['     FYI: Destination Path is ' destipath]);

% ***********************************************************
% ******************* End Pathology Choices *****************
% ***********************************************************

% read DICOMDIR
cd(sourcepath);
dcmdir = dicominfo('DICOMDIR');

% get all records in DICOMDIR
dirRecords = dcmdir.DirectoryRecordSequence;
dirRecordsFieldNames = fieldnames(dirRecords);
numOfField = size(dirRecordsFieldNames,1);

currentPatient = 0;
totalFramePerPatient = 0;

% construct a structure to hold records info from DICOMDIR
% struct dcmPatient
%           patientID
%           study
%               studyDescription
%               studyDate
%               studyTime
%               series
%                   modality
%                   seriesNumber
%                   image
%                       path
%                       numberOfFrame
for i = 1: numOfField
    currentItem = getfield(dirRecords,dirRecordsFieldNames{i});
    currentItemType = getfield(currentItem,'DirectoryRecordType');
    if strcmp(currentItemType,'PATIENT');
       if (i ~= 1)  % totalFramePerPatient has been calculated
            % display previous patient's total numbers of frames
            disp ([' Total Number of Frames: ' num2str(totalFramePerPatient)]);
       end
       
       currentPatient = currentPatient + 1;
       currentStudy = 0;
       totalFramePerPatient = 0;    %reset for each patient
       
       dcmPatient(currentPatient).patientID = currentItem.PatientID;
       disp ([num2str(currentPatient), ' Patient Name: ' currentItem.PatientName.FamilyName, ', Patient ID: ' currentItem.PatientID]);
    
    elseif strcmp(currentItemType,'STUDY')
        currentStudy = currentStudy + 1; 
        currentSeries = 0;
        
        if (isempty(currentItem.StudyDescription))
            currentItem.StudyDescription = 'No Description';
        end
        
        dcmPatient(currentPatient).study(currentStudy).studyDescription = ...
            currentItem.StudyDescription;
        dcmPatient(currentPatient).study(currentStudy).studyDate = ...
            currentItem.StudyDate;
        dcmPatient(currentPatient).study(currentStudy).studyTime = ...
            currentItem.StudyTime;
        
    elseif strcmp(currentItemType,'SERIES')
        if (currentItem.SeriesNumber ~= 0) % modality is not PR
            currentSeries = currentSeries + 1;
            currentImage = 0;
            
            dcmPatient(currentPatient).study(currentStudy).series(currentSeries).modality = ...
                currentItem.Modality;
            dcmPatient(currentPatient).study(currentStudy).series(currentSeries).seriesNumber = ...
                currentItem.SeriesNumber;
        end
       
    elseif strcmp(currentItemType,'IMAGE')
        currentImage = currentImage + 1;
       
        if (~isfield(currentItem, 'NumberOfFrames'))
            currentItem.NumberOfFrames = 1;            
        end
            
        dcmPatient(currentPatient).study(currentStudy).series(currentSeries). ...
            image(currentImage).numberOfFrames = currentItem.NumberOfFrames;
        dcmPatient(currentPatient).study(currentStudy).series(currentSeries). ...
            image(currentImage).path = currentItem.ReferencedFileID;
        
        totalFramePerPatient = totalFramePerPatient + currentItem.NumberOfFrames;      
    end
end

% display last patient's total numbers of frames
disp ([' Total Number of Frames: ' num2str(totalFramePerPatient)]);

imgInfo = dicominfo(dcmPatient(1).study(1).series(1).image(1).path);
manufacturer = imgInfo.Manufacturer;

% ask user to select one patient
% Hardcode Aug 11, 2011
% selectPatient = input('\nPlease select patient number: ');

% selectPatient = 1; % Hardcode Aug 11, 2011
disp ('Processing...');

for i=selectPatient:selectPatient   % selected patient only    
    patientID = dcmPatient(i).patientID;
    
    numOfStudy = size (dcmPatient(i).study, 2);
    for j=1:numOfStudy
        studyDescr = dcmPatient(i).study(j).studyDescription;
        studyDate = dcmPatient(i).study(j).studyDate;
        studyTime = dcmPatient(i).study(j).studyTime;
        %modality = dcmPatient(i).study(j).series(1).modality; %out of series loop, so save time
        
        numOfSeries = size (dcmPatient(i).study(j).series, 2);
        for k=1:numOfSeries
        
            modality = dcmPatient(i).study(j).series(k).modality;
            seriesNum = dcmPatient(i).study(j).series(k).seriesNumber;
          
            numOfImages = size (dcmPatient(i).study(j).series(k).image, 2);
            for m=1:numOfImages
                
                if (~isempty(dcmPatient(i).study(j).series(k).image(m)))
                    imgPath = dcmPatient(i).study(j).series(k).image(m).path;
                    imgNum = dcmPatient(i).study(j).series(k).image(m).numberOfFrames;

                    % sourcefile: sourcepath\DICOM\imagefile
                    sourcefile = [sourcepath '\' imgPath];
                    
                    % destidir:
                    % destipath\verdor\studydescr\patientID\modality\studydate_studytime\series number
                    if (diropt == 4) % special case switch for being called by parse_series
                        destidir = [destipath '\' examfolder '\' num2str(seriesNum) '\'];
                    else
                        destidir = [destipath '\' manufacturer '\' studyDescr '\' patientID '\'...
                        modality '\' studyDate '_' studyTime '\' num2str(seriesNum) '\'];
                    end % if (diropt == 4)

                    if (~exist(destidir))
                        mkdir(destidir);
                        %disp (destidir);
                    end
                    splitMultiFrame(sourcefile, destidir); % 20140211

                end
                
            end % loop of images
        end % loop of series
    end % loop of study
end %loop of patient

%disp (['Total files copied into SimpleDicom on C-Drive: ' num2str(filecount)]);
disp (' ');
disp ('Will now run "show_scan_info" on this case.  This may take awhile ... ');
%examdir = [destipath '\' series{4}{1, 1} '\' examdescr '\' cpinum '\' series{3}{1, 1}...
%            '\' examdate '_' examtime]


if (doscaninfo == 'y')
    % ********************************
    % OPTIONAL: To do a show_scan_info:
    % ********************************
    cd(destidir); % Into last series
    
    
    cd(destidir); % Into last series
    insidelastdcdir = pwd;
    cd ..; % Up one level
    outsidelastdcdir = pwd;
    lastdcdir = outsidelastdcdir;
    
    show_scan_info(1);
    % ********************************
    % ********************************
    
    if (keepmats == 'n')
        % ********************************
        % OPTIONAL: To remove image_order.mat files:
        % ********************************
        % List all series
        scanlist = dir();
        scanlist = scanlist(3:end, :);   % remove the directories . and ..

        %construct scan name (series folder) array
        scanname = {scanlist.name};  % cell array

        %convert cell array to numeric array
        numscans = size(scanname, 2);
        for i=1:numscans
            if (isdir(scanname{i}))
                thisser = scanname{i};
                cd(thisser);
                delete('*image_order.mat');
                cd ..
            end % if isdir
        end % for i
        % ********************************
        % ********************************
    end % if keepmats
end % if doscaninfo


% Create a seriesinfo.txt file.
% cd(sourcepath);
bs = findstr(sourcepath,'\'); % count backslashes.
[junk, ln2] = size(bs);
[junk, ln3] = size(sourcepath);
bss = findstr(insidelastdcdir,'\'); % count backslashes.
[junk, ln4] = size(bss);
[junk, ln5] = size(insidelastdcdir);

lastdcdir = outsidelastdcdir;
cd(lastdcdir);

% disp('Should currently be in place where show_scan_info is .. ');
% pwd
% ls
% disp('Place where seriesinfo will be stored is ')
% [sourcepath '\' sourcepath((bs(ln2)+1):ln3) '_seriesinfo.txt']


if (diropt ~= 4) % Dont recopy info only if diropt 4
    copyfile('scaninfo.txt', [sourcepath '\' sourcepath((bs(ln2)+1):ln3) '_Series_info.txt']);
end % if diropt 

cd(sourcepath);

disp('Total Script Runtime in Minutes: ');
etime(clock,tstart)/60
disp('Done');

%toc;
