function [imastor,dim1,dim2,dim3,dim4,fov1,fov2,fov3,fov4,tr,unique_te,ppd,flip,di1,nifti_params] = readdicom7(seeit, dirchoice, fpdv, subsample, slicenum, phasenum)
% NEED 'IMAGE TOOLBOX' FOR THIS MATLAB_VER13 SCRIPT.
% Function readdicom1 changes to dicom directory
% then allows one to GUI select a *.dcm file (returns 2D array,
% or to readin all *.dcm files in selected directory and return 3D array.
% If 3D array, images are sorted by escalating slice location.
% When done, return to the directory from where you came.
% Copyright 2003.  The Regents of the University of Michigan.
%
% If seeit = 1, show image(s); else dont.
% Can go to
% "http://www.gemedicalsystems.com/it_solutions/connectivity/dst.html"
% for more details on GE dicom header.

% TLChenevert,
% TLC 1/20/2003 Fix image reordering to become version "readdicom2"
% TLC 1/25/2003 Generalize for 4D data - Use "better/available" dicom tags.
%       *.InstanceNumber is in anatomics, perf, and adc GE dicom images.
%       Use it as a sort index instead of *.XXX.dcm suffice in filenames.
%       Also found *.ImagePositionPatient recognized as "better" indicator
%       of slice location.  Unfortuantely, it is a 3-vector, so for now
%       continue to use *.SliceLocation as a simple rank order of location.
%       Lastly, *.TriggerTime seems to be a reliable temporal indicator (in
%       msec), but not all series have *.TriggerTime.  Thus far, I'm unable
%       to invoke an "exist('info.TriggerTime')" command to check for
%       mulitiphase acquisitions.  So for now, let's continue to use
%       *.InstanceNumber to derive rank order.  In future, I could write a
%       dedicated readdicom for series that I know, a priori, are
%       multiphase.
% TLC 11/20/03. Generalize for file extensions w/o *.dcm(eg data from liverpool)
%               Also, in case series not in DICOM, then prompt for params
%               and input via "readimagej.m".
% TLC 12/01/03.  Minor patches to handle case when exam # = imginfo.StudyID
% is '00000'.  Then assign to 'ymmdd'.
% TLC 3/30/04.  Minor changes to return echo times for T2(*) calculation.
% TLC 4/23/04.  Fixed bug in slice ordering to insure output axial images Inferior to Superior.
% 
% HG 8/31/05    allow user to select range of images by input auguments
%               eg. readdicom5() for all, 
%                   readdicom5(3:6, 11:20) for sliced 3 to 6 of phases 11 to 20
%                   readdicom (12, 30) for slice 12 of phase 30
% TLC 9/05/05   Return TR argument as is used in rcbf scripts.
% HG 9/27/05    adjuststed to handle CT dicom files 
%
% HG 11/23/05   use function getdicom_field_num instead of imginfo.field to retrieve dicom field
%               so in the case a dicom field doesn't exist, it won't crash
% HG 12/06/05   sort ImagePositionPatient by different direction 
%                   sagittal - x
%                   coronal - y
%                   axial - z
%               if dim3*dime4 ~= numg (eg. three plane), 
%                    set dim3 = nimg, dim4 = 1
% TC 12/10/05   Switch to choose starting directory:
%                   dirchoice = 0;   % Start at usual SimpleDicom toplevel
%                   dirchoice = 1;   % Start at one level above last selection
%                   dirchoice = 2;   % Start within last selection - used to
%                   dirchoice = 3;   % Start within current directory - same as command window.
%                
%                   readin large multiphase series.
% TC 2/1/2006   Allow input from either SimpleDicom, or RadPix catchers.
% HG 7/28,2006  Each image allowed to have independent ss and rs scaling.
% TC 3/30/2007  Allow scaling to match "DV" display value on PMS scanner
%               that is for ADC and FA maps.  fpdv = 'fp' (default) returns
%               scaling that is proportional to true floating point value
%               of true MR signal.  fpdv = 'dv' returns pixel intensity
%               same as an ROI on the scanner would return.  Use "dv" for
%               reading Philips scanner-generated ADC and FA maps.
%               For fpdv = 'no' means no scaling applied.
% HG 8/14/07    add a new input parameter, subsample. if subsample == 1,
%               average pixels 2x2 to save memory
% TLC 2/20/08   Minor edit around line 228 to manually allow filenames other than "*.dcm" .
% TLC 10/28/08  Allow FOV based on "PixelSpacing" in SIEMENS dicom. Though normally, leave blotted out for speed.
% TLC 12/08/08  Generalize FOV definition in case ReconstructionDiameter is not available.
% TLC 12/15/2008 Setting to have series label 'dtes' to use series (not study) date, for multi-date exams.  
% TLC 2/10/2009 Blot-out existence check of "imginfo.StudyTime" since exist returns '0' even when it truly exists.
% TLC 11/16/2009 Add code to handle empty StudyID fields.
% TLC 1/15/2010  Add code to read ImageOrientationPatient values (direction cosines) and retain
%                them in variable "locorient" stored in image_order.mat.  Note,
%                image_order.mat files written previously by readdicom5 (and earlier) will
%                not contain variable locorient, in which case default values (0) will be assigned.
%                Also return dircosines via additional [1x6] structure element.
%                imastor.orient
% TLC 5/19/2010  Debug number of slice locs vs phases.  Retire readdicom6.
%                Denote via "MutliAngleStack Workaround" = "MASW".  Look for MASW in comments for this workaround.
% TLC 6/28/2010  Return param "di1" which simply is "dicominfo" from first
% image in imaglist, since general dicominfo can be useful for a large number of otherwise unrecorded parameters.
% TLC 7/05/2010  Some minor update for sub-dynamic pass trigger delay for Philips data to allow correct time association
%                for each individual slice in a multiphase/dynamic 2D scan.  Sub-dynamic delays NOT calculated for MRAcquisitionType='3D'.
%                FYI: the sub-dynamic delay calc'd within this script and is returned via "imastor.td" - it is NOT stored within image_order.mat.
% TLC 8/12/2010  Allow subsample = 2; for 4x4 pixel averging.
% TLC 8/24/2010  Code to utilize scale intercept and rescale intercept (sia, ria) flagged = Aug25 2010
% TLC 11/05/2010 Minor dicom file wildcard option circa line 147.
% TLC MAR 04, 2011  Allow use of ContentTime when AcquisitionTime does not
%                   exist.
% TLC Sep 14, 2011 Minor fix to check length of PixelSpacing, else show_scan_info crashes on bitmaps
% TLC Jan26_2012.  Fix for Philips dicom where "MRAcquisitionType" not defined.  Search Jan26_2012 for changes.
% TLC 20130603.  Minor bug fix on line "20130603" and in getdicom_field_str.m
% TLC 20140221.  Extend to ordering based on "ImageType" (Mag vs Real, etc), and cut-out deadwood related to "imagej".
% TLC 20140331.  Have script find suitable imgwildcard between: '*.dcm'; 'IM_*'; '0*'; '1*'; ....
% TLC 20140422.  Minor additional wildcard for HFH data 'IM0*'
% TLC 20140728:  DWI directionality. See support script "getdiffgrad.m".  Insert columns 9-12 in insortma for:
%                   insortma(9,:) = [0018,9075]: DiffusionDirectionality(0='None'; 1='Isotropic'; 2='Directional').
%                   Non-DWI sequences and DWI at b = 0 denoted '0' (None),
%                   trace DWI denoted '1' (Isotropic), directional DWI denoted '2'.
%                   [0018,9076]:
%                   DiffusionGradientDirectionSequence.Item_1.DiffusionGradientOrientation = 3vector:
%                   insortma(10,:) = cosines1 (-1 --> +1)
%                   insortma(11,:) = cosines2 (-1 --> +1)
%                   insortma(12,:) = cosines3 (-1 --> +1)
%                   insortma(13,:) = (used to be (9,:)) Holds filename index used to replace (3,:) AFTER sorting.
%                Also expand imastor structure to return DWI grad
%                directionality info, and while we're at it, image type to imastor structure too.
%                First use of this script 20140805.
% TLC 20140806:  Due to variable NumberOfAverages with b-value, need to retain image-specific NumberOfAverages in imastor structure.
%                Use insortma(13,:) to hold NumberOfAverages and insortma(14,:) for filename index.
% TLC 20140807:  Fixed bug.  Had to insert one more imgwildcard ('I*0*') for Paris data, while also catching HFH data.
%
% MD 20140820: changed some "disp" messages (for consistency with "qiba-v1" from 20140723)
% MD 20140826: commneted of "wildcard" block to use "dcmfsearch" function
% MD 20140829: commented off disp (FYIs) "Dat&", "This Series", "# images","# slices"
% MD & TLC 20140902:    Taken from QIBA SW Project
% TLC 20150106: Save nifti params in convenient structure "nifti_params"
% TLC 20150204: Put in some checks against missing fields used in
%               "nifti_params".
% TLC 20150308:  Minor check for ImageType, Rows, Columns in making nifti_params.
% TLC 20150624: No change to logic, only edited to print-out some info.
% TLC 20150627: Minor change to call "make3TRID" in creating default id name.
% TLC 20150627: Major expansion of what's stored in nifti_params. Keep full locs and orient array and more exam demographics.
% TLC 20150701: Extract for external call to "build_nifti_params.m" by both
%               readdicom7 and readdicom7_allscan.  Henceforth, expand and manage content
%               of "nifti_params" via edits to single script
%               "build_nifti_params.m".
%               Moved save('image_order',... to end of "if image_order.mat
%               exists" condition so that build_nifti_params has full
%               access to imdat information.
% TLC/DM 20150702: Much logic moved-out out of readdicom7_allscan since available here in readdicom7.
% The readdicom7_allscan script now calls readdicom7, which in turn, calls
% build_nifti_params.  We did this since we want "nifti_params" structure
% to contain all needed to create nifti and mhd files for viewing in
% 4D-Slicer (vv) and 3D-Slicer for (currently-unspecified) downstream processing.  The nifti_params structure is contained in
% "ExamDemographics.mat" (for all series in exam), as well as inside of
% "image_order.mat" (for a given series).  Keep ExamDemographics.mat around for
% future downstream processing.
% 20151105 TLC:  Insert alot of "tryGetField" calls to deal with missing
%                fields in synthetic masks from QIN DSC-Challenge datasets.
% 20151211 TLC:  Need to use "tryGetFieldFillEmpty" for TriggerTime (near line 577)

global fname dtes examID echotimes lastdcdir wordy;

warning off Images:genericDICOM;

dir0 = pwd;
%dir1 = 'C:\tlchenev_stuff\MatLab_Stuff\Misc_Data';
dir1 = dir0; % 'C:\SimpleDICOM\DicomImages_by_Vendor'; % MD 20140723
%dir1 = 'E:\SimpleDICOM\DicomImages_by_Vendor';
%dir1 = 'G:\SimpleDICOM\DicomImages_by_Vendor';
%dir1 = 'C:\Data\';

% set default values
selectall = 0;
slicestart = 1;
sliceend = 1;
phasestart = 1;
phaseend = 1;

if (isempty(wordy))
    wordy = 'y'; % Default is to print-out "Reading in image ...." which confirms new image_order.mat is being created.
end

% Set image wildcard here.  Common options are '*.dcm'; 'IM*'; '1.2*'; 'img.*' .
%imgwildcard = '*.dcm'; % default 20140331
imgwildcard = []; % 20140331 make it empty for now
%imgwildcard = '00*'; % manual override
%imgwildcard = 'img.*';
%imgwildcard = '*.IMA'; % temp manual override.

flip = [];
di1 = [];

% Regarless of input args, check if lastdcdir is defined yet - if not,
% define it.  Otherwise may get error in cd step.
ismt = isempty(lastdcdir);
if (ismt == 1)
    lastdcdir = dir1;
end % if ismt
    
if (nargin == 0)       % No input argument, readdicom5;
    seeit = 0;
    selectall = 1;
    fpdv = 'fp'; % default
    subsample = 0;
    dcdir = uigetdir(dir1, 'Pick Dicom Series Folder');
    cd(dcdir);
    lastdcdir = pwd;
elseif (nargin == 1)   % Only one input argument, readdicom5(seeit);
    selectall = 1;
    fpdv = 'fp'; % default
    subsample = 0;
    dcdir = uigetdir(dir1, 'Pick Dicom Series Folder');
    cd(dcdir);
    lastdcdir = pwd;
elseif (nargin == 2)   % Only two input args, readdicom5(seeit, dirchoice); 
    selectall = 1;
    fpdv = 'fp'; % default
    subsample = 0;
    switch dirchoice
        case 0
            dcdir = uigetdir(dir1, 'Pick Dicom Series Folder');
            cd(dcdir);
            lastdcdir = pwd;
        case 1  
            dcdir = lastdcdir;
            cd(dcdir);
            %cd ..;
            dcdir = uigetdir('', 'Pick Dicom Series Folder');
            cd(dcdir);
            lastdcdir = pwd;
        case 2
             dcdir = lastdcdir;
             cd(dcdir);
        case 3
            dir1 = pwd;
            dcdir = uigetdir(dir1, 'Pick Dicom Series Folder');
            cd(dcdir);
            lastdcdir = pwd;
        otherwise
             disp('There must be some mistake.')
    end % end switch
elseif (nargin == 3) % Only three input args, readdicom5(seeit, dirchoice, fpdv);
    selectall = 1;
    subsample = 0;
    switch dirchoice
        case 0
            dcdir = uigetdir(dir1, 'Pick Dicom Series Folder');
            cd(dcdir);
            lastdcdir = pwd;
        case 1
            dcdir = lastdcdir;
            cd(dcdir);
            %cd ..;
            dcdir = uigetdir('', 'Pick Dicom Series Folder');
            cd(dcdir);
            lastdcdir = pwd;
        case 2
             dcdir = lastdcdir;
             cd(dcdir);
        case 3
            dir1 = pwd;
            dcdir = uigetdir(dir1, 'Pick Dicom Series Folder');
            cd(dcdir);
            lastdcdir = pwd;
        otherwise
             disp('There must be some mistake.')
    end % end switch
elseif (nargin == 4) % Only four input args, readdicom5(seeit, dirchoice, fpdv, subsample);
    selectall = 1;
    switch dirchoice
        case 0
            dcdir = uigetdir(dir1, 'Pick Dicom Series Folder');
            cd(dcdir);
            lastdcdir = pwd;
        case 1
            dcdir = lastdcdir;
            cd(dcdir);
            %cd ..;
            dcdir = uigetdir('', 'Pick Dicom Series Folder');
            cd(dcdir);
            lastdcdir = pwd;
        case 2
             dcdir = lastdcdir;
             cd(dcdir);
        case 3
            dir1 = pwd;
            dcdir = uigetdir(dir1, 'Pick Dicom Series Folder');
            cd(dcdir);
            lastdcdir = pwd;
        otherwise
             disp('There must be some mistake.')
    end % end switch    
elseif ((nargin == 5) || (nargin == 6))     % Five or 6 arguments given, readdicom7(seeit, dirchoice, fpfv, subsample, slicenum, phasenum);
    slicestart = slicenum(1);
    sliceend = slicenum(end);
    phasestart = phasenum(1);
    phaseend = phasenum(end);
    switch dirchoice
        case 0
            dcdir = uigetdir(dir1, 'Pick Dicom Series Folder');
            cd(dcdir);
            lastdcdir = pwd;
        case 1
            dcdir = lastdcdir;
            cd(dcdir);
            %cd ..;
            dcdir = uigetdir('', 'Pick Dicom Series Folder');
            cd(dcdir);
            lastdcdir = pwd;
        case 2
             dcdir = lastdcdir;
             cd(dcdir);
        case 3
            dir1 = pwd;
            dcdir = uigetdir(dir1, 'Pick Dicom Series Folder');
            cd(dcdir);
            lastdcdir = pwd;
        otherwise
             disp('There must be some mistake.')
    end % end switch
else
    disp('  (FYI) Expecting: readdicom7(seeit, dirchoice, fpfv, subsample, slc1:slc2, ph1:ph2)');
end % end if nargin

thisseries = pwd;
%disp(' ');
disp(['  (FYI) Reading from: ' thisseries]);

%%% MD 20140826: commneted of "imgwildcard" block to use "dcmfsearch" function
filelist = dcmfsearch('n'); % Default.  Goes down known wildcards, then integers
%filelist = dcmfsearch('y'); % If none found, does  >>junk = dir(); and drop "." and ".."
if (isempty(filelist)) % above wild cards did not work
    disp(['ERROR: empty DICOM image file list for "' dirname '" directory ... QUITTING!']) 
    fclose(fid);% error out gracefully
    return;
end
% disp('******* Temporarily Changed filename filter1 - NEED TO UNDO IT !!! ****')
%filelist = dir('*1.3'); % Some Siemens stuff.
%filelist = dir('*IM*'); % Dicom files in this dir are stored in struct filelist.name
% *************************************************************************

nimg = length(filelist); % Also could use below: 'nimg = imginfo.Private_0025_1007(1);'

if (nimg == 0) % Perhaps not save with *.dcm file extension
    %filelist = dir('*1.2*'); % Dicom files in this dir are stored in struct filelist.name
    filelist = dir('*1.*'); % DEFAULT Dicom files in this dir are stored in struct filelist.name
    % disp('******* Temporarily Changed filename filter2 - NEED TO UNDO IT !!! ****')
    % filelist = dir('*1.3*'); % Some Siemens stuff.
    nimg = length(filelist); % Also could use below: 'nimg = imginfo.Private_0025_1007(1);'
    if (nimg == 0)
        filelist = dir('*IM*');
        nimg = length(filelist);
        if (nimg == 0)
             disp('ERROR: No image files fitting templates were found!!! Quitting ...');
             return;
        end
    end % if nimg
end % if nimg

%if (callimagej == 'n')
	% Use first in list for general image info.  Assume same for all in this series.
	di1 = dicominfo(filelist(1).name); % Save di1 in case it's called for by readdicom7 parameters
    imginfo = di1;
    if (exist('imginfo.MRAcquisitionType')) % Uncomment Jan26_2012
        mracqtype = imginfo.MRAcquisitionType; % Should be string: '2D' or '3D', some weird ones may be empty ''.
        mt_mracqtype = isempty(mracqtype);
        if(mt_mracqtype == 1)
            mracqtype = '2D'; % If mracqtype is empty, need to define something.
        end % if mt_mracqtype
    elseif (exist('imginfo.Private_2005_10a9') ) % Uncomment Jan26_2012
        mracqtype = imginfo.Private_2005_10a9; % Philips stores '2D' vs '3D' here. Jan26_2012
        mt_mracqtype = isempty(mracqtype); % Jan26_2012
        if(mt_mracqtype == 1) % Jan26_2012
            mracqtype = '2D'; % Jan26_2012
        end % if mt_mracqtype % Jan26_2012
    else % Jan26_2012
        mracqtype = 'NA';% Was Unknown Uncomment Jan25_2012
    end % if exist % Uncomment Jan26_2012
    
    % find out if the image is CT or NOT
    %isCT = strcmpi(imginfo.Modality, 'CT');
    isCT = strcmpi(getdicom_field_str(imginfo, 'Modality'), 'CT');
    
    tr = getdicom_field_num(imginfo, 'RepetitionTime');
    
    [v,d] = version;
    serDescr = getdicom_field_str(imginfo, 'SeriesDescription');
    
    if isempty(findstr(v ,'7.'))
        patName = getdicom_field_str(imginfo, 'PatientsName');
    else
        patName = getdicom_field_str(imginfo, 'PatientName');        
    end
    % disp(['This series is ' serDescr '; on Subject ' patName.FamilyName]); % TLC 20130603.
   % disp(['  (FYI) This series = "' serDescr '"; on Subject "' patName '"  ']); % MD 20140820.
    
    dim1 = getdicom_field_num(imginfo, 'Rows');
	dim2 = getdicom_field_num(imginfo, 'Columns');
    

    % Update 12/08/08. More robust definition of FOVs
    rcdiam =  getdicom_field_num(imginfo, 'ReconstructionDiameter');
    if (rcdiam == 0)
        PixelSpacing = getdicom_field_num(imginfo, 'PixelSpacing');
        fov1 = PixelSpacing(1,1)*(dim1-1);
        
        if (length(PixelSpacing) > 1)
            fov2 = PixelSpacing(2,1)*(dim2-1);
        else
            fov2 = fov1;
        end % if length Pix
    else
        fov1 = rcdiam;
        fov2 = fov1;
    end
    
    
    
%    fov1 = getdicom_field_num(imginfo, 'ReconstructionDiameter');
%    fov2 = fov1;
    
%     % Kludge setting to read SIEMENS fov:
%     PixelSpacing = getdicom_field_num(imginfo, 'PixelSpacing');
%     fov1 = PixelSpacing(1,1)*dim1;
%     fov2 = PixelSpacing(2,1)*dim2;
    
    if (subsample == 1)
        dim1 = dim1/2;
        dim2 = dim2/2;
    elseif (subsample == 2)
        dim1 = round(dim1/4);
        dim2 = round(dim2/4);
    end

    if (subsample == 1)
        pixsize = max(fov1,fov2)./double(max(dim1,dim2));
        fov1 = fov1 - pixsize;
        fov2 = fov2 - pixsize;
    elseif (subsample == 2)
        pixsize = max(fov1,fov2)./double(max(dim1,dim2));
        fov1 = fov1 - pixsize;
        fov2 = fov2 - pixsize;
    end
    
    if (isCT == 0) %MR
        %slthk = imginfo.SpacingBetweenSlices; %  Note this incluces slice thickness plus gap.
        slthk = getdicom_field_num(imginfo, 'SpacingBetweenSlices');
    else %CT
        %slthk = imginfo.SliceThickness 
        slthk = getdicom_field_num(imginfo, 'SliceThickness');
    end
    
    % if manufacturer is Philips, get scaling slope and prepause delay
    %manufacturer = imginfo.Manufacturer;
    manufacturer = getdicom_field_str(imginfo, 'Manufacturer');
    if isempty(findstr(manufacturer ,'Philips'))
        isPhilips = 0;
        ss = 1;
        rs = 1;
        ppd = 0;
    else % Philips
        isPhilips = 1;
        %[ss, rs] = getscale(imginfo);
        [ss, rs, si, ri] = getscale(imginfo);
        ppd = getppd(imginfo);
    end
    
    % TLC 12/01/2003 Patch to fix HFH anonymized data format that sets
    % exno=00000.
    
    % *****************************************************
    % Start TLC fixes for Empty Study cells Nov16, 2009 ***
    % *****************************************************
    excell = imginfo.StudyID;
    if (isempty(excell))
        imginfo.StudyID = '9999';
    end
    excell = imginfo.StudyTime;
    if (isempty(excell))
        imginfo.StudyTime = '9999';
    end    
    
    examnumber = str2num(imginfo.StudyID);
    % End  TLC fixes Nov16, 2009 ***
    
    if (examnumber == 0) % As HFH anonymizer does
        imginfo.StudyID = imginfo.StudyDate(4:8);
    end
    
	examID = ['e' imginfo.StudyID]; % A char array
    
	% Example: dtes = 'd20020101010101e4138s1'
	%dtes = ['d' yy mm dd hhmm(1:2) 'e' num2str(exno) 's' num2str(serno)];
    
	%if (exist('imginfo.StudyTime') == 0) % Have to check this for HFH data
    %    imginfo.StudyTime = '00'; 
    %end
    
    
    % 12/15/2008:
    % May choose to use SERIES DATE instead of default EXAM DATE here for studies that span multiple scandates:
    %dtes = ['d' imginfo.SeriesDate '00e' imginfo.StudyID 's' num2str(imginfo.SeriesNumber)];
    %imginfo.StudyTime(1:2);
    
	%dtes = ['d' imginfo.StudyDate imginfo.StudyTime(1:2) 'e' imginfo.StudyID 's' num2str(imginfo.SeriesNumber)];  % 20140707.
    sdate = tryGetField(imginfo,'SeriesDate','00000000'); % 20150204 returns "0000000000" if no field present
    fndum4 = tryGetField(imginfo,'SeriesTime','0000'); % 20150204 returns "0000" if no field present
    stime = fndum4(1:4); % 20150204
    % dtes = ['d' imginfo.SeriesDate imginfo.StudyTime(1:2) 'e' imginfo.StudyID 's' num2str(imginfo.SeriesNumber)]; % 20140707. 
    dtes = ['d' sdate stime(1:2) 'e' imginfo.StudyID 's' num2str(imginfo.SeriesNumber)]; % 20150204
	%disp(['  (FYI) Date & ExamID Info = "' dtes '"  ']); % MD: 20140820
	%disp(['  (FYI) # of Images in Series = "' num2str(nimg) '"  ']);

    % ***********************************************************************
    % ***********************************************************************
    % ************ image_order.mat does NOT exist - make one ****************
    % ***********************************************************************
    % ***********************************************************************
    
    if (exist('image_order.mat', 'file') == 0) % Following implemented if image_order.mat does NOT exist.
        posxord = zeros(1,nimg); % Will be used to discriminate image position order, x direction.
        posyord = zeros(1,nimg); % Will be used to discriminate image position order, y direction.
        poszord = zeros(1,nimg); % Will be used to discriminate image position order, z direction.
        % MutliAngleStack Workaround = "MASW"
        ao = [-800 -900 -1000]; % MASW. Arb origin skewed. Makes output order of simple 3-plane to be (Ax,Sag,Cor)
        posrord = zeros(1,nimg); % MASW.  This will be the radial distance from an arbitrary origin 'ao'
        
        instord = zeros(1,nimg); % Will be used to discriminate acquisition order.
        te = zeros(1,nimg); % echo time array
        %flip = zeros(1,nimg); % TLC 7/8/08
        
        if (isPhilips == 1)
            ssa = zeros(1, nimg); % scale slope array
            rsa = zeros(1, nimg); % rescale slope array
            sia = zeros(1, nimg); % scale intercept array
            ria = zeros(1, nimg); % rescale intercept array
        end % isPhilips
        
        TagTD = isfield(imginfo,'TriggerTime');
          
        for ii = 1:nimg
            %fprintf('.'); % Shows progress 'dots' as each image read-in.
            if (wordy == 'y')
                disp([' Reading in image = ' num2str(ii)]);
            end
            imginfo = dicominfo(filelist(ii).name);
            
            posexist = isfield(imginfo, 'ImagePositionPatient');
            if (posexist == 0)
                posxord(ii) = 0;  %default
                posyord(ii) = 0;  %default
                poszord(ii) = 0;  %default
            else
                posxord(ii) = imginfo.ImagePositionPatient(1,1);
                posyord(ii) = imginfo.ImagePositionPatient(2,1);
                poszord(ii) = imginfo.ImagePositionPatient(3,1);
                % MutliAngleStack Workaround = "MASW"
                posrord(ii) = sqrt(((ao(1)-posxord(ii))^2) + ((ao(2)-posyord(ii))^2) + ((ao(3)-poszord(ii))^2)); % MASW
                
            end % if posexist
            
            %instord(ii) = imginfo.InstanceNumber;
            instord(ii) = getdicom_field_num(imginfo, 'InstanceNumber');
            
            fileord{ii} = filelist(ii).name;    %store file names in the order of input
            
            %te(ii) = getdicom_field_num(imginfo, 'EchoTime');
            te(ii) = tryGetFieldFillEmpty(imginfo,'EchoTime',0); % 20151105 TLC
            flip(ii) = getdicom_field_num(imginfo, 'FlipAngle'); % TLC 7/8/08
            
            if (isPhilips == 1)
                %[ssa(ii), rsa(ii)] = getscale(imginfo);
                [ssa(ii), rsa(ii), sia(ii), ria(ii)] = getscale(imginfo);
            end % isPhilips
            
            % construct insortma matrix (nimg x 8)
            % 1st column is a place holder for ImagePositionPatient
            insortma(ii,1) = 0;
            
            % 2nd column, trigger time / instance number
            if (TagTD == 1)
                
%                 % Start kludge ...
%               %filelist(ii).name
%               if (isfield(ii,'TriggerTime') == 0)
%                   insortma(ii,2)= 0;
%               else
%                   insortma(ii,2)= imginfo.TriggerTime;
%               end % if (isfield(ii,'TriggerTime')
%               % end kludge.  Unblot next line when nokludge
              %insortma(ii,2) = imginfo.TriggerTime; % 20151211
              insortma(ii,2) = tryGetFieldFillEmpty(imginfo, 'TriggerTime', 0); % 20151211
              
            else
              insortma(ii,2) = instord(ii);
            end % if TagTD
           
            % 20140221.  3rd column, track the index of image data
            % insortma(ii,3) = ii; % effectively this is the alphanumeric
            % order of filenames.  Lets put in a dummy "1" in whole 3rd column
            % so it does not drive sorting since filename is arbitrary.
            % HOWEVER, to maintain backward compatibility with existing
            % image_order.mat files on disk, the 3rd column needs to remain
            % filename index, so AFTER "sortrows", restore 3rd column with
            % (appropriatly-sorted) filename index.
            insortma(ii,3) = 1; % For now 3rd column is nonfunctional. 
            
            % 4th column, echo time
            insortma(ii,4) = te(ii);

            % MAR 04, 2011.  Sometimes, 'AcquisitionTime' does not exist,
            % but 'ContentTime' does.
            acqexist = isfield(imginfo, 'AcquisitionTime');
            if (acqexist == 0)
                acqexist1 = isfield(imginfo, 'ContentTime');
                if (acqexist1 == 0)
                    insortma(ii,5) = 0;  %default
                else
                    insortma(ii,5) = 1000*(str2num(imginfo.ContentTime(1:2))*3600 + str2num(imginfo.ContentTime(3:4))*60 + str2num(imginfo.ContentTime(5:end)));
                end % if acqexist1
            else
                 insortma(ii,5) = 1000*(str2num(imginfo.AcquisitionTime(1:2))*3600 + str2num(imginfo.AcquisitionTime(3:4))*60 + str2num(imginfo.AcquisitionTime(5:end)));
            end % if acqexist
                
            % 6th column, b value
            insortma(ii,6) = getbvalue(imginfo); % 20140707
            %insortma(ii,6) = 0; % 20140707
            
            % 7th column, ppd
            insortma(ii,7) = getppd(imginfo);
            
            
            % 8th column for image type, eg mag, adc, real, imagingary etc. An 20140221 addition. 
            insortma(ii,8) = getimagetype(imginfo); % this returns a number 0,1,2,... based on imagetype
            
%             % 9th column:  Eventually, 9th (& 3rd column) will contain 
%             insortma(ii,9) = ii; % Keep track of filename index "ii" to confirm 3rd col restored correctly.
            
            % 20140728 Start ... Get Diffusion Gradient Info:
            [diffdir, diffdircosines] = getdiffgrad(imginfo);
            insortma(ii,9) = diffdir; % 0=None; 1=Isotropic; 2=Directional
            insortma(ii,10) = diffdircosines(1);
            insortma(ii,11) = diffdircosines(2);
            insortma(ii,12) = diffdircosines(3);
            % 20140728 ... End
            % 20140806 13th column:
            navexist = isfield(imginfo, 'NumberOfAverages');
            if (navexist == 0) % No number of averaged
                insortma(ii,13) = 0; % Number of averages = 0 will be code for 'not found'
            else
                insortma(ii,13) = imginfo.NumberOfAverages;
            end % if navexist
            % 20140806 End.
            
            % 14th column:
            insortma(ii,14) = ii; % Keep track of filename index "ii" to confirm 3rd col restored correctly.
            

            % location matrix (3 x nimg)
            if (posexist == 0)
                locpos(1,ii) = 0;  %default
                locpos(2,ii) = 0;  %default
                locpos(3,ii) = 0;  %default
                locorient(1,ii) = 1; %axial default;
                locorient(2,ii) = 0; %axial default;
                locorient(3,ii) = 0; %axial default;
                locorient(4,ii) = 0; %axial default;
                locorient(5,ii) = 1; %axial default;
                locorient(6,ii) = 0; %axial default;
            else
                locpos(1,ii) = imginfo.ImagePositionPatient(1,1);  
                locpos(2,ii) = imginfo.ImagePositionPatient(2,1);  
                locpos(3,ii) = imginfo.ImagePositionPatient(3,1);
                locorient(1,ii) = imginfo.ImageOrientationPatient(1,1);
                locorient(2,ii) = imginfo.ImageOrientationPatient(2,1);
                locorient(3,ii) = imginfo.ImageOrientationPatient(3,1);
                locorient(4,ii) = imginfo.ImageOrientationPatient(4,1);
                locorient(5,ii) = imginfo.ImageOrientationPatient(5,1);
                locorient(6,ii) = imginfo.ImageOrientationPatient(6,1);
            end % if posexist

       end % for ii
       
       %fprintf('\r'); % print a new line
       
       %retrieve unique elements of te array
       unique_te = unique(te);
       
       %get position difference of each direction
       posxdiff = max(posxord) - min(posxord);
       posydiff = max(posyord) - min(posyord);
       poszdiff = max(poszord) - min(poszord);
       
       posrdiff = max(posrord) - min(posrord) - 2; % MASW.  "-2" favors posxyzdiff over posrdiff, since it avoids a tie.
       
       %get the biggest difference among x, y, z directions
       % posdiffmax = max([posxdiff, posydiff, poszdiff]); blotted MASW
       posdiffmax = max([posxdiff, posydiff, poszdiff, posrdiff]); % MASW
       
       insortma(:, 1) = []; %delete 1st column, it was just a place holder
       
       if (posxdiff == posdiffmax) % sort position by x direction
           insortma = [posxord', insortma];
           [sortpos isortpos] = sort(posxord);
           %disp (['x direction']);
       elseif (posydiff == posdiffmax) % sort position by y direction
           insortma = [posyord', insortma];
           [sortpos isortpos] = sort(posyord);
           %disp (['y direction']);
       elseif (poszdiff == posdiffmax) % sort position by z direction
           insortma = [poszord', insortma];
           [sortpos isortpos] = sort(poszord);
           %disp (['z direction']);
       elseif (posrdiff == posdiffmax) % MASW
           insortma = [posrord', insortma]; % MASW
           [sortpos isortpos] = sort(posrord); % MASW
           %disp (['r direction']); % MASW         
       end % if posxdiff
        
       % Now calculate number of unique image locations:
       dsortpos = abs(sortpos(2:nimg) - sortpos(1:(nimg-1))); % A derivative-like array.
       nlocs = sum(ones(1,(nimg-1)) & floor(dsortpos/0.00001)) + 1;   %Normalize non-zero derivatives; sum; add 1; ans = nlocs.
       dim3 = nlocs; % number of image locations
       
       % *************** MANUAL OVERIDE NUMBER OF LOCS ******************
       % Need to manually fix "2-stack series with 2-overlapping slices"
       % Delete image_order.mat, set overide values, re-run readdicom5.
       % dim3 = 24;
       % nlocs = dim3;
       % *************** END MANUAL OVERIDE *****************************     
       
       dim4 = ceil(nimg/nlocs); % number of phases
          
       if (dim3*dim4 ~= nimg) %eg. three plane data % MD: 20140820 chnaged warning to "disp" message
           disp(' WARNING: dim3(# slices) * dim4(# phases) not equal to (# images)... Resetting dim3 and dim4');
           dim3 = nimg;
           dim4 = 1;
       end % if dim3*dim4
       
       %output number of slices and number of phases
       %disp(['   (FYI): # of Slices = "' num2str(dim3) '"  ']);
       %disp(['  (FYI): # of Phases = ' num2str(dim4)]);
        
       % Sort 3rd & 4th dimension order by image position, instancenumber (
       % or triggertime when existed ), serially by columns in insortma.
       % Highest sort priority is 1st column, 2nd column is next sort
       % priority.
       % Last column is "ii", filename index - need this to track and relate to
       % fileord and ssa, rsa, sia, ria.
       % Use "isortedre" to resort other indexed lists, fileord, locpos, locorient, ssa, etc..
       [sortedre,isortedre] = sortrows(insortma);
       
       % ***********************************************
       % ***********************************************
       %% Debug point - compare pre- and post-sort arrays
       % save('prepostsort','insortma','sortedre');
       % ***********************************************
       % ***********************************************
       
       %sortedre = sortrows(insortma);
       % Now that sorting is done, restore filename index in column3 to
       % maintain backward compatibility with prior image_order.mat files.
%        for ii = 1:nimg
%            sortedre(ii,3) = isortedre(ii);
%        end % for ii
       % The above loop works fine, but easier to use this:
       % Now that sorting is done (where filename index was the weakest sort criterion),restore fn index in column 3 used for read-in order.
       % sortedre(:,3) = sortedre(:,13);
       sortedre(:,3) = sortedre(:,14); % 20140806 Update to 14th column.
       
       %save order info to file image_order.mat
       % ****** "sortedre" is now sorted in logical ascending
       % order based on slice location (column1), triggertime (col2), NOT by
       % filename index (col3), ... YES by imagetype (col8), then YES by filename index
       % (col9). NOTE!!!! fileord, locpos, locorient, ssa, rsa, sia, ria ARE NOT
       % sorted in-sync with sortedre!!!  They WILL BE sync'd with sortedre
       % upon reading-in the actual images thanks to indexing fileord etc based on
       % "sortedre(:,3)".  See this in nested loops "for i" and "for j" below.
       % Avoid changing things in these nested loops below else risk loss
       % of backward compatibility.
  
       %set image selection borders
       if (selectall == 1)
           sliceend = dim3;
           phaseend = dim4;
       end
       
       % build a struct matrix(ixj) to hold return values
       tic % tic toc a pair commands to measure elipsed time
       
       % Setup "acqord" and "subpassfly" which is used by Philips for 2D slice interleaving. Define acqord all as 1's for 3D.
       nslice = dim3;
       if (mracqtype == '3D')
            acqord = ones(1,nslice);
            subpassdly = acqord; % ie all slices at the same time
       else
            nstep = round(sqrt(nslice));
            %nstep = 2; % 20130321 TEMP kludge to force odd-slice even-slice interleaving
            [indx, acqord] = ph_interleave(nslice,nstep);
            [indx, subpassdly] = sort(acqord,2);
       end % if mracqtype
     
       for i = slicestart:sliceend
            for j = phasestart:phaseend
                if (((i-1)*dim4 + j) <= nimg)  % there are enough images, matrix dimensions not exceeded
                    % retrieve image data based on sorted orders
                    xx = dicomread(fileord{sortedre((i-1)*dim4+j,3)});

                    if (isPhilips == 0) %non-Philips data, no scaling
                        x = double(xx);
                    else
                        ss = ssa(sortedre((i-1)*dim4+j,3));
                        rs = rsa(sortedre((i-1)*dim4+j,3));
                        si = sia(sortedre((i-1)*dim4+j,3));  % Aug25 2010
                        ri = ria(sortedre((i-1)*dim4+j,3));  % Aug25 2010
                        %disp(['ScaleSlope = ' num2str(ss) '; RescaleSlope = ' num2str(rs)]);
                        if (fpdv == 'dv')
                            % x = double(xx) * rs);
                            x = (double(xx) * rs) + ri; % Aug25 2010
                        elseif (fpdv == 'fp')
                            % x = (1e-2)*double(xx/(ss*rs)); %For image scaling per Philips.
                            %  x = double(xx)/ss; % PREFERRED For image scaling.
                            
                            %disp(['ss = ' num2str(ss) ';  si = ' num2str(si)]);
                            
                              x = (double(xx) - si)/ss;  % Aug25 2010
                        else % then must mean 'no' scaling
                            % x = double(xx);
                             x = double(xx) - si;  % Aug25 2010
                        end % if fpdv
                        
                    end % if isPhilips
                    clear xx;
                    
                    if (subsample == 1)
                       x = reduceImg2(x); %average pixels 2x2, save memory
                    elseif (subsample == 2)
                       x = reduceImg4(x); %average pixels 4x4, save memory
                    end % if subsample

                    imastor(i-slicestart+1,j-phasestart+1).idata = x;
                    clear x;
                    
                    % Rev 7/06/2010 make "if" more specific to TagTD and Philips vs nonPhilips.
                    % OLD WAY: if (TagTD == 1) % trigger delay is available
                    if ( (TagTD == 1) && (isPhilips == 0) ) % trigger delay is available and not Philips data
                       imastor(i-slicestart+1,j-phasestart+1).td =  sortedre((i-1)*dim4+j,2);
                    else % Philips data
                      if (dim4 == 1) % only one phase, set trigger delay as 0
                          imastor(i-slicestart+1,j-phasestart+1).td = 0;
                      else %multple phases

                          % Calc sub-Dynamic pass acquisition delay for each slice in 2D multiphase Philips images
                          if (j == dim4) %last phase
                            %imastor(i-slicestart+1,j-phasestart+1).td = sortedre(j,5) + (acqord(i)-1)*(sortedre(j,5) - sortedre(j-1,5))/dim3 - sortedre(1,5);
                            imastor(i-slicestart+1,j-phasestart+1).td = sortedre(j,5) + (((subpassdly(i)-1)*(sortedre(j,5) - sortedre(j-1,5)))/dim3) - sortedre(1,5);
                          else 
                            %imastor(i-slicestart+1,j-phasestart+1).td = sortedre(j,5) + (acqord(i)-1)*(sortedre(j+1,5) - sortedre(j,5))/dim3 - sortedre(1,5);
                            imastor(i-slicestart+1,j-phasestart+1).td = sortedre(j,5) + (((subpassdly(i)-1)*(sortedre(j+1,5) - sortedre(j,5)))/dim3) - sortedre(1,5);
                          end % j == dim4
                      end %dim4 == 1
                    end %TagTD == 1  

                    imastor(i-slicestart+1,j-phasestart+1).loc =  locpos(:,sortedre((i-1)*dim4+j,3));
                    imastor(i-slicestart+1,j-phasestart+1).echotime = sortedre((i-1)*dim4+j,4);
                    imastor(i-slicestart+1,j-phasestart+1).bvalue = sortedre((i-1)*dim4+j,6);
                    imastor(i-slicestart+1,j-phasestart+1).ppd = sortedre((i-1)*dim4+j,7);
                    imastor(i-slicestart+1,j-phasestart+1).orient  = locorient(:,sortedre((i-1)*dim4+j,3));
                    
                    % 20140728.  Add diffusion gradient directionality, and image type in structure:
                    imastor(i-slicestart+1,j-phasestart+1).diffdir = [sortedre((i-1)*dim4+j,9); sortedre((i-1)*dim4+j,10); sortedre((i-1)*dim4+j,11); sortedre((i-1)*dim4+j,12)];
                    imastor(i-slicestart+1,j-phasestart+1).imgtyp = sortedre((i-1)*dim4+j,8);
                    % 20140806.  Retain NumberOfAverages in structure:
                    imastor(i-slicestart+1,j-phasestart+1).naverages = sortedre((i-1)*dim4+j,13);
                    
                end %if (((i-1)*dim4 + j) <= nimg)
            end %for j
       end  %for i
       
       fov3 = (dim3 - 1)*slthk;
       fov4 = dim4; % Dont really know this parameter yet.
       
       %show images
       numphase = phaseend-phasestart+1;
       numslice = sliceend-slicestart+1;
       
       if seeit == 1
            for i =1:numphase
                for j=1:numslice
%                     junk2 = imastor(j, i).loc;
%                     disp('XYZ values:');
%                     junk2'
%                     junk2 = imastor(j, i).orient;
%                     disp('DIRCosine values:');
%                     junk2'
                    img(imastor(j, i).idata);
                    pause;
                end
            end
       end % if seeit
       
      % Moved build_nifti_params and save('image_order.mat' ...) to here since
      % now imastor is fully created, thus easy to get info for nifti_params from it:
      [nslc nphase] = size(imastor); % This reliably returns [nslc nphase] but for some fields, getsafield transposes dependent on nslc=1 vs nslc >1.
%       % GEOMETRY STUFF:
%       junk = getsafield(imastor,'loc');
%       nip.loc = junk(:,:,1); % Assume 1st vector plane reps locs for all.
%       junk = getsafield(imastor,'orient');
%       nip.orient = junk(:,:,1); % Assume 1st vector plane reps orient for all.
      
      % GEOMETRY and NON-GEOMETRY STUFF, but possibly 4th-dimension:
      if (nphase > 1)
          
          if (nslc == 1) % [1 nphase] = size(imastor);
              junk = (getsafield(imastor,'td'))';       % junk = ([nphase 1])' = [1 nphase]
              nip.td = junk;                            % [1 nphase]
              junk = (getsafield(imastor,'echotime'))'; % junk = ([nphase 1])' = [1 nphase]
              nip.echotime = junk;                      % [1 nphase]
              junk = (getsafield(imastor,'ppd'))';      % junk = ([nphase 1])' = [1 nphase]
              nip.ppd = junk;                           % [1 nphase]
              junk = (getsafield(imastor,'imgtyp'))';   % junk = ([nphase 1])' = [1 nphase]
              nip.imgtyp = junk;                        % [1 nphase]
              junk = (getsafield(imastor,'naverages'))';% junk = ([nphase 1])' = [1 nphase]
              nip.naverages = junk;                     % [1 nphase]
              junk = (getsafield(imastor,'bvalue'))';   % junk = ([nphase 1])' = [1 nphase]
              nip.bvalue = junk;                        % [1 nphase]
              junk = getsafield(imastor,'diffdir');     % junk = [4 nphase]
              nip.diffdir = junk;                       % [4 nphase]
              junk = getsafield(imastor,'loc');         % junk = [3 nphase]
              nip.loc = junk(:,1);                      % [3 1] only has 1 slice, same loc all phases
              junk = getsafield(imastor,'orient');      % junk = [6 nphase]
              nip.orient = junk(:,1);                   % [6 1] only has 1 slice, same orient all phases
          else % then multi-slice, multi-phase
              junk = getsafield(imastor,'td');          % junk =  [nslc nphase]
              nip.td = junk;                            % [nslc nphase] keep all slice info
              junk = getsafield(imastor,'echotime');    % junk =  [nslc nphase]
              nip.echotime = junk(1,:);                 % [1 nphase] only need 1 slice
              junk = getsafield(imastor,'ppd');         % junk =  [nslc nphase]
              nip.ppd = junk;                           % [nslc nphase] keep all slice info
              junk = getsafield(imastor,'imgtyp');      % junk =  [nslc nphase]
              nip.imgtyp = junk(1,:);                   % [1 nphase] only need 1 slice
              junk = getsafield(imastor,'naverages');   % junk =  [nslc nphase]
              nip.naverages = junk(1,:);                % [1 nphase] only need 1 slice
              junk = getsafield(imastor,'bvalue');      % junk =  [nslc nphase]
              nip.bvalue = junk(1,:);                   % [1 nphase] only need 1 slice
              junk = getsafield(imastor,'diffdir');     % junk = [4 nslc nphase]
              nip.diffdir = squeeze(junk(:,1,:));       % [4 nphase] only need 1 slice
              junk = getsafield(imastor,'loc');         % junk = [3 nslc nphase]
              nip.loc = squeeze(junk(:,:,1));           % [3 nslc] keep all slice loc info, same loc all phases 
              junk = getsafield(imastor,'orient');      % junk = [6 nslc nphase]
              nip.orient = squeeze(junk(:,:,1));        % [6 nslc] keep all slice orient info, same orient all phases
          end % if nslc == 1
          
      else % then only single phase, but assume (for now) multi-slice
          
          junk = (getsafield(imastor,'td'))';          % junk =  ([1 nslc])' = [nslc 1]
          nip.td = junk;                               % [nslc 1] keep all slice info
          junk = getsafield(imastor,'echotime');       % junk =  [1 nslc]
          nip.echotime = junk(1,1);                    % [1 1] only need 1 slice
          junk = (getsafield(imastor,'ppd'))';         % junk =  ([1 nslc])' = [nslc 1]
          nip.ppd = junk;                              % [nslc 1] keep all slice info    
          junk = (getsafield(imastor,'imgtyp'))';      % junk =  ([1 nslc])' = [nslc 1]
          nip.imgtyp = junk(1,1);                      % [1 1] only need 1 slice
          junk = getsafield(imastor,'naverages');      % junk =  [1 nslc]
          nip.naverages = junk(1,1);                   % [1 1] only need 1 slice
          junk = getsafield(imastor,'bvalue');         % junk =  [1 nslc]
          nip.bvalue = junk(1,1);                      % [1 1] only need 1 slice
          junk = getsafield(imastor,'diffdir');        % junk =  [4 nslc]
          nip.diffdir = junk(:,1);                     % [4 1] only need 1 slice
          junk = getsafield(imastor,'loc');            % junk =  [3 nslc]
          nip.loc = junk;                              % [3 nslc] keep all slice info
          junk = getsafield(imastor,'orient');         % junk =  [6 nslc]
          nip.orient = junk;                           % [6 nslc] keep all slice info
          
      end % if nphase > 1
      
      % nifti_params = build_nifti_paramsWIP(imginfo,locpos,locorient,sortedre,fov1,fov2,dim3,dim4);
      nifti_params = build_nifti_params(imginfo,nip,fov1,fov2,dim3,dim4);

      if (isPhilips == 0)
          save('image_order.mat', 'sortedre', 'dim3', 'dim4', 'fileord', 'locpos', 'locorient','TagTD', 'tr', 'unique_te', 'ppd', 'nifti_params'); % As of 20150105 same nifti_params too.
      else
          save('image_order.mat', 'sortedre', 'dim3', 'dim4', 'fileord', 'locpos', 'locorient','TagTD', 'tr', 'unique_te', 'ppd', 'nifti_params','ssa', 'rsa', 'sia', 'ria'); % As of 20150105 same nifti_params too.   
      end % isPhilips 
       
     
       
    % ***********************************************************
    % ***********************************************************
    % *********** image_order.mat does exist - use it ***********
    % ***********************************************************
    % ***********************************************************
    
    else % Following implemented if image_order DOES exist

        %retrieve order info from image_order.mat
        load ('image_order')
        
        % Need to account for image_order.mat files created before readdicom6 where locorient does not exist.
        existorient = exist('locorient','var');
        if (existorient == 0) % locorient does not exist, set to zero & warn user.
            disp(' ');
            disp(' WARNING:  Image orientation Info not saved in this image_order.mat');
            disp('           Assigned as straight axial ... ');
            disp('           If orientation needed, delete image_order.mat from series folder and (re)run.');
            disp(' ');
            locorient = zeros(6,(dim3*dim4));
            locorient(1,:) = ones(1,(dim3*dim4));
            locorient(5,:) = ones(1,(dim3*dim4));
        end % if existorient
        
        %output number of slices and number of phases
        %disp(['  (FYI) # of Slices = "' num2str(dim3) '"  ']);
        % disp(['  (FYI): # of Phases = ' num2str(dim4)]);
       
        %set image selection borders
        if (selectall == 1)
           sliceend = dim3;
           phaseend = dim4;
        end % if selectall
        
        % build a struct matrix(ixj) to hold return values
        tic % tic toc a pair commands to measure elipsed time
        
        % Setup "acqord" and "subpassfly" which is used by Philips for 2D slice interleaving. Define acqord all as 1's for 3D.
        nslice = dim3;
        if (mracqtype == '3D')
            acqord = ones(1,nslice);
            subpassdly = acqord; % ie all slices at the same time
        else
            nstep = round(sqrt(nslice));
            [indx, acqord] = ph_interleave(nslice,nstep);
            [indx, subpassdly] = sort(acqord,2);
        end % if mracqtype
        

        for i = slicestart:sliceend
            for j = phasestart:phaseend
                %disp(['Slice ' num2str(i) '; Phase ' num2str(j)]);
                if (((i-1)*dim4 + j) <= nimg)  % there are enough images, matrix dimensions not exceeded
                   % retrieve image data based on sorted orders
                   
%                    % Mar 15, 2011 DTI debug stuff start *****
%                    dii = dicominfo(fileord{sortedre((i-1)*dim4+j,3)});
%                    dgo = dii.DiffusionGradientOrientation;
%                    [-dgo(2) dgo(1) -dgo(3)]/0.7071 % Need this to make print-out look like DTI_medium_overplus.txt
%%                                     |  0   -1    0  |     dgo(1)
%%                    This implies:    |  1    0    0  |  x  dgo(2)
%%                                     |  0    0   -1  |     dgo(3)
%%                    % Mar 15, 2011 DTi debug stuff end *****
                   
%                    sortedre((i-1)*dim4+j,3) % 20150624
%                    fileord{sortedre((i-1)*dim4+j,3)} % 20150624
                   xx = dicomread(fileord{sortedre((i-1)*dim4+j,3)});
                   
%                    % TEMP Hack to debug variable Naverages vs b-value:
%                    dixx = dicominfo(fileord{sortedre((i-1)*dim4+j,3)});
%                    [fileord{sortedre((i-1)*dim4+j,3)} '    ' num2str(dixx.NumberOfAverages)]

                   if (isPhilips == 0) %non-Philips data, no scaling
                        x = double(xx);
                   else    % Philips data, need scaling
                        if ( exist('ssa', 'var') == 1)
                            ss = ssa(sortedre((i-1)*dim4+j,3));
                        end
                        
                        if ( exist('rsa', 'var') == 1)   
                            rs = rsa(sortedre((i-1)*dim4+j,3));
                        end
                        
                        if ( exist('sia', 'var') == 1) % Aug25 2010 
                            si = sia(sortedre((i-1)*dim4+j,3));
                        end
                        
                        if ( exist('ria', 'var') == 1) % Aug25 2010 
                            ri = ria(sortedre((i-1)*dim4+j,3));
                        end
                        
                        
                        if (fpdv == 'dv')
                            % x = double(xx) * rs;
                            x = (double(xx) * rs) + ri; % Aug25 2010
                        elseif (fpdv == 'fp')
                            % x = (1e-2)*double(xx/(ss*rs)); %For image scaling per Philips.
                            % x = double(xx)/ss; % PREFERRED For image scaling.
                             x = (double(xx) - si)/ss;  % Aug25 2010
                        else % must be 'no'
                            % x = double(xx); % No scaling
                             x = double(xx) - si;  % Aug25 2010
                        end % if fpdv
                        
                   end % if isPhilips
                   clear xx;

                   if (subsample == 1)
                       x = reduceImg2(x); %average pixels 2x2, save memory
                   elseif (subsample == 2)
                       x = reduceImg4(x); %average pixels 4x4, save memory
                   end % if subsample
                   
                   imastor(i-slicestart+1,j-phasestart+1).idata = x;
                   clear x;
                   
                   
                   % Rev 7/06/2010 make "if" more specific to TagTD and Philips vs nonPhilips.
                   % OLD WAY: if (TagTD == 1) % trigger delay is available
                   if ( (TagTD == 1) && (isPhilips == 0) ) % trigger delay is available and not Philips data
                       imastor(i-slicestart+1,j-phasestart+1).td =  sortedre((i-1)*dim4+j,2);
                   else % Philips data
                      if (dim4 == 1) % only one phase, set trigger delay as 0
                          imastor(i-slicestart+1,j-phasestart+1).td = 0;
                      else %multple phases

                          % Calc sub-Dynamic pass acquisition delay for each slice in 2D multiphase Philips images
                          if (j == dim4) %last phase
                            %imastor(i-slicestart+1,j-phasestart+1).td = sortedre(j,5) + (((acqord(i)-1)*(sortedre(j,5) - sortedre(j-1,5)))/dim3) - sortedre(1,5);
                            imastor(i-slicestart+1,j-phasestart+1).td = sortedre(j,5) + (((subpassdly(i)-1)*(sortedre(j,5) - sortedre(j-1,5)))/dim3) - sortedre(1,5);
                          else 
                            %imastor(i-slicestart+1,j-phasestart+1).td = sortedre(j,5) + (((acqord(i)-1)*(sortedre(j+1,5) - sortedre(j,5)))/dim3) - sortedre(1,5);
                            imastor(i-slicestart+1,j-phasestart+1).td = sortedre(j,5) + (((subpassdly(i)-1)*(sortedre(j+1,5) - sortedre(j,5)))/dim3) - sortedre(1,5);
                          end % j == dim4
                          
                      end %dim4 == 1
                   end %TagTD == 1  

                   imastor(i-slicestart+1,j-phasestart+1).loc = locpos(:,sortedre((i-1)*dim4+j,3));
                   imastor(i-slicestart+1,j-phasestart+1).echotime = sortedre((i-1)*dim4+j,4);
                   imastor(i-slicestart+1,j-phasestart+1).bvalue = sortedre((i-1)*dim4+j,6);
                   imastor(i-slicestart+1,j-phasestart+1).ppd = sortedre((i-1)*dim4+j,7);
                   imastor(i-slicestart+1,j-phasestart+1).orient = locorient(:,sortedre((i-1)*dim4+j,3));
                   
                   % 20140728.  Add diffusion gradient directionality:
                   imastor(i-slicestart+1,j-phasestart+1).diffdir = [sortedre((i-1)*dim4+j,9); sortedre((i-1)*dim4+j,10); sortedre((i-1)*dim4+j,11); sortedre((i-1)*dim4+j,12)];
                   imastor(i-slicestart+1,j-phasestart+1).imgtyp = sortedre((i-1)*dim4+j,8);
                   % 20140806.  Add NumberOfAverages:
                   imastor(i-slicestart+1,j-phasestart+1).naverages = sortedre((i-1)*dim4+j,13);
                    
                end % if (((i-1)*dim4 + j) <= nimg)  % there are enough images, matrix dimensions not exceeded
            end % for j=phasestart:phaseend
        end % for i=slicestart:sliceend
                
        fov3 = (dim3 - 1)*slthk;
        fov4 = dim4; % Dont really know this parameter yet.
        
        %show images
        numphase = phaseend-phasestart+1;
        numslice = sliceend-slicestart+1;
       
        if seeit == 1
            for i =1:numphase
                for j=1:numslice
                    img(imastor(j, i).idata);
                    pause;
                end
            end
        end % if seeit
    
    end  % if image_order.mat exists  
    % ************************************************
    % ******** done with image_order.mat logic *******
    % ************************************************

t = toc;

cd(dir0); % Return to starting directory.

%end % End function call
