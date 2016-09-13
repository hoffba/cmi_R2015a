function sortscans = show_scan_info(browz)
% display info of all scans of selected exam as following example
%   Series Number = 14
%   Patient Name = ...
%   Series Description = DWI B=100
%   Scan Sequence = EP\SE
%   Repetition Time = 0
%   Unique Echo Time = 102.9
%   Prepause Delay = NaN
%   Number of Slices = 24
%   Number of Phases = 5
%   Imaging Frequency = 638901280
% 11/16/05   HG

% browz = 0; User has to browse down from top level dicom directory;
% browz = 1; Browse from last selected location
% TLC Jan 2010.  Based on "readdicom6_allscan" to handle orientation volume info.
% TLC Apr 2010: Minor work-around to handle Carlo Pierpaoli series named "002", "003"
%               ... as opposed to "2","3". See transient fix around line 71.
%               After you are done, UNDO the fix and put back to normal.
% TLC May 2010: Minor change to readdicom7_allscan around line 93.
% TLC Aug 11, 2011: Include In-plane FOV and Image Dimesionality in txt file.
% TLC 20130604 in anticipation of Matlab 2013
% TLC 20130903 Also list series times and durations.
% TLC 20140221 Test call to readdicom7_allscantlc().  Seems to work,
% therefore retired then replaced as "readdicom7_allscan"
% TLC 20140331  Generalize so this script can be run directly on
%               TRAID-created exam folders, without having to go into and
%               out of KPacs.  Search for " 20140331 " for changes.
%               Create structures for SiteScanner demographics and
%               ExamSeries parameters so that we can sort according to
%               SeriesNumber as opposed to FolderNumber (unrelated to
%               series number).  Once sorted, the scaninfo.txt will be
%               created.
% TLC 20140502  Create wildwildcard workaround option. DONT FORGET TO
% CHANGE wildwildcard to 'y' IN BOTH THIS SCRIPT AND readdicom7_allscan.
% TLC 20140807:  Fixed bug.  Had to insert one more imgwildcard ('I*0*') for Paris-Kpacs export data, while also catching HFH data.
%
% MD 20140820: add "qiba_split_series" call
% MD & TLC 20140902:    Taken from QIBA SW Project
% TLC 20141223: Bug fix for SECONDARY captures/scanned docs that do not have legit SeriesTime.  magtime no longer used so blot deadwood.
% TLC 20150106: Save a few more params used in nifti creation in ExamDemographics.
% TLC/DM 20150703: Call to new readdicom7_allscan that calls, readdicom7.
% TLC 20150706: Get scanner demographics from middle-ish series (1st-few and last-few are often screenshots).
% TLC 20150716: Still need greater immunity to missing standard dicom
%               fields like "Institution" etc.  Use alot more "tryGetField.m".
% TLC 20151211: Change unknown SeriesTime default from 'UNK' to '000000', to handle screenshots that do not have sertime.
% TLC 20160520: Save to JUNK.mat for debugging only - no substantative changes.
% TLC 20160520: Changes to "FillEmpty" SeriesNumber (Siemens).
% MD 20160523:  extract additional "group/dna" params into ExamSeries:
%               pixbw,freqmx, phasemx, parallel, pifactforinplane
% MD 20160523: added a call to "get_vendorCode" and "get_pifactor" (for parallel imaging info)
% MD 20160523: changed "ExamSeries" structure assignement to "tempstruct"
% TLC 20160523: Switch to " if (isequal({some_string},{'target_string'})) " to check equality of variable length srings. 
% TLC 20160524: Alternative to 20160523 change. Instead use MatLab's suggestion: "if( strcmp(str1,str2) )".

global fname dtes examID lastdcdir; % This is so I can label images with appropriate filenames.

%select exam of interest
if (nargin == 0)       % No input argument,
    browz=0; % Assign browse starting point to default if browz is empty
    disp('  *** Select EXAM FOLDER containing all the series.    ***');
    disp('  (FYI) Do NOT go so deep as to select a SERIES FOLDER!  ');
    %examdir = uigetdir('E:\SimpleDICOM\DicomImages_by_Vendor', 'Pick Exam of Interest');
    examdir = uigetdir('C:\data', 'Pick Exam of Interest');
else
    examdir = lastdcdir;
end

disp('Starting show_scan_info ...');
tstart = clock;

% CHANGE wildwildcard to 'y' IN BOTH THIS SCRIPT AND readdicom7_allscan.
wildwildcard = 'n'; % 20140502. def = 'n'.  'y' means image filenames found by wide-open "*" wildcard.  Be sure to manually delete all extraneous files (infos, jpgs, mats,...).
%examdir
cd(examdir);

scanlist = dir();
scanlist = scanlist(3:end, :);   % remove the directories "." and ".." (assumes "Windows")

chkdir = getsafield(scanlist, 'isdir')'; % MD 20140820

if (length(chkdir) == length(chkdir(chkdir==0)) && (length(chkdir(chkdir==0))>3)) % no series folders MD 20140820
    % try to split series into folders
    disp(' ')
    disp('  (FYI) Selected DICOM data is NOT SPLIT into series: SPLITTING ............ ');
    examdir = split_series(examdir); % update "examdir" to the one with split series
    cd(examdir); % MD 20140820: enter series directory before cataloging
    %pwd
    if (exist('scaninfo.txt', 'file') || exist('ExamDemographics.mat', 'file')) % MD: 20140722
        disp(' ');
        %disp('  (FYI) Found existing  catalog files for selected exam: OVERWRITING ............');
        %MD 20140829: NOTE existing catalog is not overwritten
        disp('  (FYI) Found existing  catalog files for selected exam: USING ............'); 
%         disp(' ');
%         disp('  (FYI) To RE-CATALOG: remove existing "scaninfo.txt" and "ExamDemographics.mat" ')
%        disp(['  from current DICOM DIR = "' examdir '"'])
        % return; % MD 20140826: will use existing catalog files
        %%% otherwise, for default re-cataloging could use the code below:
        % delete('scaninfo.txt', 'file'); delete('ExamDemographics.mat','file');
    end
    scanlist = dir(); % MD 20140820: update scan-list for series
    scanlist = scanlist(3:end, :);   % remove the directories "." and ".." (assumes "Windows")
    chkdir = getsafield(scanlist, 'isdir'); % MD 20140820
end


[numscans junk] = size(scanlist); % 20140331
midseries = ceil(numscans/2);

ser1img1 = []; % 20140331 This will retain dicominfo of 1st series, 1st image

%%% MD 20140826: to use existing "scaninfo" do check here:
%%% if (~exist('scaninfo.txt', 'file'))
fid = fopen('scaninfo.txt', 'wt'); % will "overwrite" existing catalog without above "ifexist" check
fprintf(fid, '%s\n\n', examdir);

% 20141223 Blot next line
t = []; % Time in magnet in minutes

folderindex = 0; % Count folders 1,2,3...
%loop thru each scan and display scan info
serieslist = [];
disp(['  (FYI) Cataloging detected ' num2str(length(chkdir(chkdir==1))) ' scan series ............'])
%for i = 1:(numscans-4) % Use when last few series folders are empty or other non-standard dicom format.
for i = 1:numscans
    
    if ( scanlist(i).isdir == 1) % 20140331
        folderindex = folderindex + 1;

        % Next line is default:
        % dirname = num2str(sortedscans(i)); % 20140331
        dirname = scanlist(i).name; % 20140331

        % But in case you need to use the text name because someone named series folders "002", "003", ...
        % dirname = char(scanname(i));

        cd(dirname);

        % *********************************************************************
        % *********************************************************************
        % dcmlist=dir('*.dcm'); % 20140331 blot. DEFAULT Dicom files in this dir are stored in struct filelist.name
        %disp('******* Temporarily Changed filename filter3 - NEED TO UNDO IT !!! ****')
        %dcmlist = dir('00*'); % Some Siemens stuff.
        % *********************************************************************
        % *********************************************************************
        
        % Start 20140331 ..
        % Filter down to dicom files only, using high probability filters first.
        dcmlist = dcmfsearch(wildwildcard); % MD 20140826: replaced the code-block below with a function

        if (isempty(dcmlist)) % above wild cards did not work
            disp(['ERROR: empty DICOM image file list for "' dirname '" directory ... QUITTING!']) 
            fclose(fid);% error out gracefully
            return;
        end
        % Save dicominfo from 1st series 1st file, since for Siemens at
        % least, the first series has orignal MR images whereas last series
        % tends to be derived/screenshots/not-full dicom stuff.
        % if ( isempty(ser1img1) )
        if (i == midseries) % Using middle series instead to avoid misleading scanner demographics from non-MRI scanner devices
            % disp(['MidSeries DirName is ' dirname]); % 20160520 debug
            fname = dcmlist(1).name;
            ser1img1 = dicominfo(fname);
            % Put site & scanner demographics on top o  scaninfo.txt
            tempval = tryGetField(ser1img1,'InstitutionName','UNK');
            %fprintf(fid, 'InstitutionName = %s\n', getdicom_field_str(ser1img1, 'InstitutionName'));
            fprintf(fid, 'InstitutionName = %s\n', tempval);
            SiteScanDemographics.Institution = tempval;
            
            tempval = tryGetField(ser1img1,'Manufacturer','UNK');
            %fprintf(fid, 'Manufacturer = %s\n', getdicom_field_str(ser1img1, 'Manufacturer'));
            fprintf(fid, 'Manufacturer = %s\n', tempval);
            SiteScanDemographics.Manufacturer = tempval;
            
            tempval = tryGetField(ser1img1,'ManufacturerModelName','UNK');
            % fprintf(fid, 'Model = %s\n', getdicom_field_str(ser1img1, 'ManufacturerModelName'));
            fprintf(fid, 'Model = %s\n', tempval);
            SiteScanDemographics.Model = tempval;
            
            tempval = tryGetField(ser1img1,'MagneticFieldStrength',0);
            %fprintf(fid, 'FieldStrength (T) = %4.2f\n', getdicom_field_num(ser1img1,'MagneticFieldStrength'));
            fprintf(fid, 'FieldStrength (T) = %4.2f\n', tempval);
            SiteScanDemographics.FieldStrength = tempval;
            
            tempval = tryGetField(ser1img1,'StationName','UNK');
            %fprintf(fid, 'StationName = %s\n', getdicom_field_str(ser1img1, 'StationName'));
            fprintf(fid, 'StationName = %s\n', tempval);
            SiteScanDemographics.StationName = tempval;

            tempval = tryGetField(ser1img1,'DeviceSerialNumber','UNK');
            %fprintf(fid, 'DeviceSerialNumber = %s\n', getdicom_field_str(ser1img1, 'DeviceSerialNumber'));
            fprintf(fid, 'DeviceSerialNumber = %s\n', tempval);
            SiteScanDemographics.SerialNumber = tempval;
            
            tempval = tryGetField(ser1img1,'SoftwareVersion','UNK');
            %fprintf(fid, 'SoftwareVersion = %s\n', getdicom_field_str(ser1img1, 'SoftwareVersion'));
            fprintf(fid, 'SoftwareVersion = %s\n', tempval);
            SiteScanDemographics.SoftwareVersion = tempval;
            
            tempval = tryGetField(ser1img1,'StudyDescription','UNK');
            %fprintf(fid, 'StudyDescription = %s\n', getdicom_field_str(ser1img1, 'StudyDescription'));
            fprintf(fid, 'StudyDescription = %s\n', tempval);
            SiteScanDemographics.StudyDescription = tempval;
            
            tempval = tryGetField(ser1img1,'StudyDate','UNK');
            %fprintf(fid, 'StudyDate = %s\n', getdicom_field_str(ser1img1, 'StudyDate'));
            fprintf(fid, 'StudyDate = %s\n', tempval);
            SiteScanDemographics.StudyDate = tempval;
            
            tempval = tryGetField(ser1img1,'StudyTime','UNK');
            %fprintf(fid, 'StudyTime = %s\n', getdicom_field_str(ser1img1, 'StudyTime'));
            fprintf(fid, 'StudyTime = %s\n', tempval);
            SiteScanDemographics.StudyTime = tempval;
            
            tempval = tryGetField(ser1img1,'PatientName','UNK');
            SiteScanDemographics.PatientName = tempval;
            
            fprintf(fid, 'FYI: Nifti & MHD Params Saved \n\n'); % 20150702, note that MHD info saved too
            
            % For "check_filesize" script, we need patient name listed in scaninfo.txt after each series number
            % So blot next line and insert PatientName in loops below.
            % fprintf(fid, 'Patient Family Name = %s\n\n',getdicom_field_str(ser1img1, 'PatientName')); 
            % Build SiteScanner demographic structure:
            % SiteScanDemographics.Institution = getdicom_field_str(ser1img1, 'InstitutionName');
            % SiteScanDemographics.Manufacturer = getdicom_field_str(ser1img1, 'Manufacturer');
            % SiteScanDemographics.Model = getdicom_field_str(ser1img1, 'ManufacturerModelName');
            % SiteScanDemographics.FieldStrength = getdicom_field_str(ser1img1, 'MagneticFieldStrength');
            % SiteScanDemographics.StationName = getdicom_field_str(ser1img1, 'StationName');
            % SiteScanDemographics.SerialNumber = getdicom_field_str(ser1img1, 'DeviceSerialNumber');
            % SiteScanDemographics.SoftwareVersion = getdicom_field_str(ser1img1, 'SoftwareVersion');
            % SiteScanDemographics.StudyDescription = getdicom_field_str(ser1img1, 'StudyDescription');
            % SiteScanDemographics.StudyDate = getdicom_field_str(ser1img1, 'StudyDate');
            % SiteScanDemographics.StudyTime = getdicom_field_str(ser1img1, 'StudyTime');
            % SiteScanDemographics.PatientName = getdicom_field_str(ser1img1, 'PatientName');
        end % if i == midseries
        % ... End 20140331
            
        fname = dcmlist(1).name;
        imginfo=dicominfo(fname);
        
        tempval = tryGetField(imginfo,'SeriesDescription','UNK');
        %seriesdescr = getdicom_field_str(imginfo, 'SeriesDescription');
        seriesdescr = tempval;

        % tempval = tryGetField(imginfo,'SeriesNumber',0);        % 20160520
        tempval = tryGetFieldFillEmpty(imginfo,'SeriesNumber',0); % 20160520
        
        %seriesnumber = getdicom_field_num(imginfo, 'SeriesNumber'); %20140331
        seriesnumber = tempval;
        
        tempval = tryGetField(imginfo,'ScanningSequence','UNK');
        %scanseq = getdicom_field_str(imginfo, 'ScanningSequence');
        scanseq = tempval;
        
        %tempval = tryGetField(imginfo,'PatientName');
        %patientname = getdicom_field_str(imginfo, 'PatientName');
        %patientname = tempval;
        %patfamname = patientname.FamilyName;
        %patfamname = patientname; % 20130604 TLC
        %patfamname = tempval;
        
        tempval = tryGetField(imginfo,'ImagingFrequency',0);
        %imgfreq = getdicom_field_num(imginfo, 'ImagingFrequency');
        imgfreq = tempval;

        % 20130903 TLC.  Show info related to scan time & durations.  Start ...
        
        %tempval = tryGetField(imginfo,'SeriesTime','UNK'); % 20151211
        tempval = tryGetField(imginfo,'SeriesTime','000000'); % 20151211
        %sertime = getdicom_field_str(imginfo, 'SeriesTime');
        sertime = tempval;
 
% 20141223 Start BLOT ... Bug fix for SECONDARY captures & scanned-docs that do not have legit SeriesTime.
%        sertsecs = 60*60*str2num(sertime(1:2)) + 60*str2num(sertime(3:4)) + str2num(sertime(5:6));
        
        
%         if (isempty(sertsecs))
%             sertsecs = 0;
%         else
%             t = [t sertsecs/60];
%         end % if esempty
%         
%         % ... End 20141223
%         
%         magtime = (sertsecs/60) - t(1);
% 20141223 ... End BLOT Bug fix for SECONDARY captures & scanned-docs that do not have legit SeriesTime

        if (isfield(imginfo, 'AcquisitionDuration'))
            tempval = tryGetField(imginfo,'AcquisitionDuration',0);
            %serduration = getdicom_field_num(imginfo, 'AcquisitionDuration');
            serduration = tempval;
        elseif (isfield(imginfo, 'FrameAcquisitionDuration'))
            tempval = tryGetField(imginfo,'FrameAcquisitionDuration',0);
            %serdurationmsecs = getdicom_field_num(imginfo, 'FrameAcquisitionDuration'); % in msec
            serdurationmsecs = tempval;
            serduration = serdurationmsecs/1000;
        else
            serduration = 0; % Not found
        end
        % serduration = round(serdurationflt);
        % ... End 20130903 TLC.


        %get info from readdicom5_allscan()    
        %[imastor,dim1,dim2,numslice,numphase,fov1,fov2,fov3,fov4,tr,unique_te,ppd] = readdicom5_allscan();
        %[dim1,dim2,numslice,numphase,fov1,fov2,tr,unique_te,ppd] = readdicom5_allscan();
        %[dim1,dim2,numslice,numphase,fov1,fov2,tr,unique_te,ppd,flip] = readdicom6_allscan();
        dims = [1 1 1 1];
        fov = [0 0];
        %[dim1,dim2,numslice,numphase,fov1,fov2,tr,unique_te,ppd,flip] = readdicom7_allscan();
        %disp('calling readdicom7_allscantlc ...');
        %[dim1,dim2,numslice,numphase,fov1,fov2,tr,unique_te,ppd,flip] = readdicom7_allscantlc();
        %[dim1,dim2,numslice,numphase,fov1,fov2,tr,unique_te,ppd,flip] = readdicom7_allscan();
        %[dim1,dim2,numslice,numphase,fov1,fov2,tr,unique_te,ppd] = readdicom7_allscan(); % 20140331
        %[dim1,dim2,numslice,numphase,fov1,fov2,tr,unique_te,ppd] = readdicom7_allscan(); % 20150106
        %[dim1,dim2,numslice,numphase,fov1,fov2,tr,unique_te,ppd] = readdicom7_allscan_WIP(); % 20150204
        %[dim1,dim2,numslice,numphase,fov1,fov2,tr,unique_te,ppd] = readdicom7_allscan(); % 20150204
        %[dim1,dim2,numslice,numphase,fov1,fov2,tr,unique_te,ppd] = readdicom7_allscanWIP(); % 20150627
        %[dim1,dim2,numslice,numphase,fov1,fov2,tr,unique_te,ppd] = readdicom7_allscan(); % 20150627
        %[dim1,dim2,numslice,numphase,fov1,fov2,tr,unique_te,ppd] = readdicom7_allscanBakTLCWIP(); % 20150627
        [dim1,dim2,numslice,numphase,fov1,fov2,tr,unique_te,ppd] = readdicom7_allscan(); % Much leaner readdicom7_allscan 20150703
        
        fov = [fov1 fov2];
        dims = [dim1 dim2 numslice numphase];

        %fprintf('\n\n'); % print a new line

        % write scan info into file scaninfo.txt
        % fprintf(fid, 'Series Number = %u\n', sortedscans(i)); % 20140331
%         fprintf(fid, 'Series Folder = %s\n', dirname); % 20140331
%         fprintf(fid, 'Series Number = %u\n', seriesnumber); % 20140331
%         fprintf(fid, 'Patient Family Name = %s\n', patfamname);
%         fprintf(fid, 'Series Description = %s\n', seriesdescr);
%         fprintf(fid, 'Scan Sequence = %s\n', scanseq);
%         fprintf(fid, 'Repetition Time = %g\n', tr);
%         fprintf(fid, 'Unique Echo Time = %g\n', unique_te);
%         fprintf(fid, 'Prepause Delay = %g\n', ppd);
%         fprintf(fid, 'Number of Slices %u\n', numslice);
%         fprintf(fid, 'Number of Phases = %u\n', numphase);
%         fprintf(fid, 'Image Dimensionality = %u %u %u %u\n', dims);
%         fprintf(fid, 'In-plane FOV (mm) = %6.1f %6.1f\n', fov);
%         %fprintf(fid, 'Imaging Frequency = %g\n', imgfreq);                     % 20140331 TLC
%         fprintf(fid, 'Series Start in Military Time = %s\n', sertime(1:6));     % 20130903 TLC
%         fprintf(fid, 'Relative Start Time of Series (min) = %f\n', magtime);    % 20130903 TLC
%         fprintf(fid, 'Series Duration (secs) = %f\n\n', serduration);           % 20130903 TLC

% MD 20160523: extract additional "group/dna" params into ExamSeries:
    pxbw = tryGetField(imginfo,'PixelBandwidth',0);
    acqmx = tryGetField(imginfo,'AcquisitionMatrix',[0 0 0 0]);
    ippedir = tryGetField(imginfo,'InPlanePhaseEncodingDirection', 'UNK');
    vendor = get_vendorCode(tryGetField(imginfo,'Manufacturer', 'UNK'));
    %vendor  % debug
    [prl pfinpln] = get_pifactor(imginfo, vendor);
    %if (ippedir == 'ROW') % eg [0;508;319;0] = 0\508\319\0 % 20160523
    %if (isequal({ippedir},{'ROW'})) % 20160523
    if (strcmp(ippedir,'ROW')) % 20160524
        frqmx = acqmx(2);  % New for ACRIN6698
        phsmx = acqmx(3); % New for ACRIN6698
    %elseif (ippedir == 'COL') % then is 'COL' or 'UNK' % eg [160;0;0;180] = 160\0\0\180 % 20160523
    %elseif (isequal({ippedir},{'COL'})) % then is 'COL' or 'UNK' % eg [160;0;0;180] = 160\0\0\180 % 20160523
    elseif (strcmp(ippedir,'COL')) % then is 'COL' or 'UNK' % eg [160;0;0;180] = 160\0\0\180 % 20160524
        frqmx = acqmx(1);
        phsmx = acqmx(4);
    else % UNK
       frqmx = 0;
       phsmx = 0;
    end % if Group.ippedir
    tmpstruct(folderindex).pixbw = pxbw;
    tmpstruct(folderindex).freqmx = frqmx;
    tmpstruct(folderindex).phasemx = phsmx;
    tmpstruct(folderindex).parallel = prl;
    tmpstruct(folderindex).pifactorinplane = pfinpln;
% MD 20160523: END extract additional 
        serieslist = [serieslist;seriesnumber];
        tmpstruct(folderindex).Folder = dirname;
        tmpstruct(folderindex).SeriesNumber = seriesnumber;
        tmpstruct(folderindex).SeriesDescription = seriesdescr;
        tmpstruct(folderindex).Sequence = scanseq;
        tmpstruct(folderindex).TR = tr;
        tmpstruct(folderindex).TEs = unique_te;
        tmpstruct(folderindex).PrepauseDelay = ppd;
        tmpstruct(folderindex).NSlice = numslice;
        tmpstruct(folderindex).NPhase = numphase;
        tmpstruct(folderindex).ImageDimensions = dims;
        tmpstruct(folderindex).InPlaneFOV = fov;
        tmpstruct(folderindex).SerStartTime = sertime(1:6);
        tmpstruct(folderindex).SerDurationSecs = serduration;
        
        % 20150106 Start ... Get nitfi params from image_order.mat
        load image_order;
        tmpstruct(folderindex).NiftiParams = nifti_params;
        tmpstruct(folderindex).SpacingBetweenSlices = nifti_params.slicectr2ctr;
        % ... End 20150106

        %go back to the exam directory
        cd ..;
    end % if (scanlist(i).isdir == 1) % 20140331


end % for numscans

[sortseries, isort] = sort(serieslist);
% save('C:\temp\JUNK','serieslist','sortseries','isort','tmpstruct','folderindex'); % Debug save 20160520

for iser = 1:folderindex
    % MD 20160523: changed structure assignment to single-step
    ExamSeries(iser) = tmpstruct(isort(iser));
%     ExamSeries(iser).SeriesNumber = tmpstruct(isort(iser)).SeriesNumber;
%     ExamSeries(iser).Folder = tmpstruct(isort(iser)).Folder;
%     ExamSeries(iser).SeriesDescription = tmpstruct(isort(iser)).SeriesDescription;
%     ExamSeries(iser).Sequence = tmpstruct(isort(iser)).Sequence;
%     ExamSeries(iser).TR = tmpstruct(isort(iser)).TR;
%     ExamSeries(iser).TEs = tmpstruct(isort(iser)).TEs;
%     ExamSeries(iser).PrepauseDelay = tmpstruct(isort(iser)).PrepauseDelay;
%     ExamSeries(iser).NSlice = tmpstruct(isort(iser)).NSlice;
%     ExamSeries(iser).NPhase = tmpstruct(isort(iser)).NPhase;
%     ExamSeries(iser).ImageDimensions = tmpstruct(isort(iser)).ImageDimensions;
%     ExamSeries(iser).InPlaneFOV = tmpstruct(isort(iser)).InPlaneFOV;
%     ExamSeries(iser).SerStartTime = tmpstruct(isort(iser)).SerStartTime;
%     ExamSeries(iser).SerDurationSecs = tmpstruct(isort(iser)).SerDurationSecs;
%     
%     % 20150106 Start ... Save nitfi params from image_order.mat
%     ExamSeries(iser).NiftiParams = tmpstruct(isort(iser)).NiftiParams; % As of 20150702, NiftiParams will have info for retrospective MHD creation (4D and Demographics too)
%     ExamSeries(iser).SpacingBetweenSlices = tmpstruct(isort(iser)).SpacingBetweenSlices;
    % ... End 20150106 
       
    % Finish scaninfo.txt file, but sorted by increasing series number:
    fprintf(fid, 'Series Number = %u\n', ExamSeries(iser).SeriesNumber);
    fprintf(fid, 'Series Folder = %s\n', ExamSeries(iser).Folder);
    % Redundant, but need Patient Name repeated for each series to keep "check_filesize" correct.
    
    tempval = tryGetField(ser1img1,'PatientName','UNK');
    %fprintf(fid, 'Patient Family Name = %s\n', getdicom_field_str(ser1img1, 'PatientName'));
    fprintf(fid, 'Patient Family Name = %s\n', tempval);
    fprintf(fid, 'Series Description = %s\n', ExamSeries(iser).SeriesDescription);
    fprintf(fid, 'Scan Sequence = %s\n', ExamSeries(iser).Sequence);
    fprintf(fid, 'Repetition Time = %g\n', ExamSeries(iser).TR);
    fprintf(fid, 'Unique Echo Time = %g\n', ExamSeries(iser).TEs);
    fprintf(fid, 'Prepause Delay = %g\n', ExamSeries(iser).PrepauseDelay);
    fprintf(fid, 'Number of Slices %u\n', ExamSeries(iser).NSlice);
    fprintf(fid, 'Slice Ctr-Ctr Spacing (mm) = %6.2f\n', ExamSeries(iser).SpacingBetweenSlices);
    fprintf(fid, 'Number of Phases = %u\n', ExamSeries(iser).NPhase);
    fprintf(fid, 'Image Dimensionality = %u %u %u %u\n', ExamSeries(iser).ImageDimensions);
    fprintf(fid, 'In-plane FOV (mm) = %6.1f %6.1f\n', ExamSeries(iser).InPlaneFOV);
    fprintf(fid, 'Series Start in Military Time = %s\n', ExamSeries(iser).SerStartTime);
    fprintf(fid, 'Series Duration (secs) = %f\n\n', ExamSeries(iser).SerDurationSecs);
end % for iser

%disp('lastdcdir leaving show_scan_info is ');
% lastdcdir
% close file

fclose(fid);
save('ExamDemographics','examdir','SiteScanDemographics','ExamSeries');
%%% end    % if no "scaninfo.txt" 
disp(['  Total show_scan_info runtime = ' num2str(etime(clock,tstart)/60) ' minutes.']);
%etime(clock,tstart)/60

return;