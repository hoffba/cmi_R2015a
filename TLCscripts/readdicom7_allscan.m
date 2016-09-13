function [dim1,dim2,dim3,dim4,fov1,fov2,tr,unique_te,ppd,flip] = readdicom7_allscan()
%function [imastor,dim1,dim2,dim3,dim4,fov1,fov2,fov3,fov4,tr,unique_te,ppd] = readdicom6_allscan();
% modified version of readdicom5()
% called from show_scan_info()
% no input augument needed
% no image shown, no need to GUI select series
%
% HG 7/28,2006  Each image allowed to have independent ss and rs scaling.
% HG 8/10,2007  remove return varibles imastor, fov3 and fov4 to save time
%               and memory.
% TLC 7/21/2008 Drop #secs in 'imginfo.AcquisitionTime', since may not be present in GE
% series.
% TLC 1/15/2010  Add code to read ImageOrientationPatient values (direction cosines) and retain
%                them in variable "locorient" stored in image_order.mat.  Note,
%                image_order.mat files written previously by readdicom5_allscan (and earlier) will
%                not contain variable locorient, in which case default values (0) will be assigned.
%                Also return dircosines via additional [1x6] structure element.
%                imastor.orient.  There is a corollary "readdicom6.m"
% TLC 5/19/2010  Debug number of slice locs vs phases.  Retire readdicom6_allscan.
%                Denote via "MutliAngleStack Workaround" = "MASW".  Look for MASW in comments for this workaround.
% TLC 8/25/2010  Minor change to utilize (re)scale intercepts if they exist flag = Aug25 2010.
% TLC MAR 04, 2011  Allow use of ContentTime when AcquisitionTime does not
%                   exist.
% DIM July 22, 2011:  Fix logic in check for TagTD near line 187.
% TLC Aug 11, 2011: Make fov1 & fov2 definitions same as readdicom7 - more robust. Search by "Aug 11, 2011".
% TLC Sep 14, 2011: SAVE sia & ria arrays that were debugged on Aug 25, 2010.  Debugged, but were not saved!
%                   Also fixed check on PixelSpacing for bitmaps.
% TLC 20140221.  Parallel new readdicom7. Extend to ordering based on "ImageType" (Mag vs Real, etc), and cut-out deadwood related to "imagej".
% TLC 20140331.  Have script find suitable dicom imagefile wildcard
% TLC 20140422.  Minor additional wildcard for HFH data 'IM0*'
% TLC 20140502.  wildwildcard workaround option.
% TLC 20140728.  DWI directionality. See support script "getdiffgrad.m".  Insert columns 9-12 in insortma for:
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
%                directionality info, and while we're at it, add image type to imastor structure too.
%                Replaced previous version for first use on 20140805.
% TLC 20140806.  Retain image-specific NumberOfAverages in imastor structure.  insortma(13,:) now holds NumberOfAverages and insortma(14,:) for filename index.
% 
% TLC 20140807:  Fixed bug.  Had to insert one more imgwildcard ('I*0*') for Paris data, while also catching HFH data.
% MD 20140820: chnaged warning/error "disp" messages
% MD & TLC 20140902:    Taken from QIBA SW Project
% TLC 20150106: Save nifti params in convenient structure "nifti_params"
% TLC 20150204: Put in some checks against missing fields used in
%               "nifti_params".
% TLC 20150308:  Minor check for ImageType, Rows, Columns in making nifti_params.
% TLC 20150627:  More info into nifti_params.
% TLC 20150629: Extract for external call to "build_nifti_params.m" by both
%               readdicom7 and readdicom7_allscan.  Henceforth, expand and manage content
%               of "nifti_params" via edits to single script
%               "build_nifti_params.m".
% TLC/DM 20150703: Much logic moved-out since available in readdicom7.
% This script now calls readdicom7, which in turn, calls
% build_nifti_params.  We did this since we want "nifti_params" structure
% to contain all needed to create nifti and mhd files for viewing in
% 4D-Slicer (vv) and 3D-Slicer for (currently-unspecified) downstream processing.  The nifti_params structure is contained in
% "ExamDemographics.mat" (for all series in exam), as well as inside of
% "image_order.mat" (for a given series).  Keep ExamDemographics.mat around for
% future downstream processing. 

global lastdcdir wordy
warning off Images:genericDICOM;

% CHANGE wildwildcard to 'y' IN BOTH THIS SCRIPT AND readdicom7_allscan.
wildwildcard = 'n'; % 20140502. def = 'n'.  'y' means image filenames found by wide-open "*" wildcard.  Be sure to manually delete all extraneous files (infos, jpgs, mats,...).

% set default values
slicestart = 1;
sliceend = 1;
phasestart = 1;
phaseend = 1;
seeit = 0;
selectall = 1;

filelist = dcmfsearch(wildwildcard); 
if (isempty(filelist)) % above wild cards did not work
    disp(['ERROR: empty DICOM image file list for "' dirname '" directory ... QUITTING!']) 
    fclose(fid);% error out gracefully
    return;
end

nimg = length(filelist); % Also could use below: 'nimg = imginfo.Private_0025_1007(1);'

if (nimg == 0) % Perhaps not save with *.dcm file extension
    filelist = dir('*1.*'); % Dicom files in this dir are stored in struct filelist.name
    nimg = length(filelist); % Also could use below: 'nimg = imginfo.Private_0025_1007(1);'
    if (nimg == 0)
          disp('ERROR: No image files fitting templates were found!!! Quitting ...');
          return;
    end % if nimg
end % if nimg

imginfo = dicominfo(filelist(1).name);

isCT = strcmpi(getdicom_field_str(imginfo, 'Modality'), 'CT');

tr = getdicom_field_num(imginfo, 'RepetitionTime');
dim1 = getdicom_field_num(imginfo, 'Rows');
dim2 = getdicom_field_num(imginfo, 'Columns');

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
    ppd = NaN;
else % Philips
    isPhilips = 1;
    %[ss, rs] = getscale(imginfo);
    [ss, rs, si, ri] = getscale(imginfo); % Aug25 2010
    ppd = getppd(imginfo);
end

% TLC 12/01/2003 Patch to fix HFH anonymized data format that sets
% exno=00000.
examnumber = str2num(imginfo.StudyID);
if (examnumber == 0) % As HFH anonymizer does
    imginfo.StudyID = imginfo.StudyDate(4:8);
end

examID = ['e' imginfo.StudyID]; % A char array

% Example: dtes = 'd20020101010101e4138s1'
if (exist('imginfo.StudyTime') == 0) % Have to check this for HFH data
    imginfo.StudyTime = '00'; 
end

dtes = ['d' imginfo.StudyDate imginfo.StudyTime(1:2) 'e' imginfo.StudyID 's' num2str(imginfo.SeriesNumber)];

if (exist('image_order.mat', 'file') == 0) 
    locdir = pwd;
    lastdcdir = locdir;
    browz = 2; % choose "current" directory
    wordy = 'n'; % pass this to "readdicom7" to avoid image-by-image "print-out" when readdicom7 called from readdicom7_allscan.
    %%% create "image_order" and "nifit_params" for this series
    [imdat,dim1,dim2,dim3,dim4,fov1,fov2,fov3,fov4,tr,unique_te,ppd,flip,dinf, nifti_params] = readdicom7(0,browz);
    cd(locdir);
    wordy = 'y'; % turn wordy back on for stand-alone readdicom7 calls.
else
    locdir = pwd;
    [pathstr,thisname,ext] = fileparts(locdir);
    lastdcdir = locdir;
    disp(['    FYI: Series folder ' thisname ' already has an image_order.mat.  Using it ... ']);
    cd(locdir);
    wordy = 'y'; % turn wordy back on for stand-alone readdicom7 calls.
end % if exist

load ('image_order'); % Some output arguments were stored in "image_order.mat" so need to load it.
