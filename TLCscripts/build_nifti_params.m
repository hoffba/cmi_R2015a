function [ nifti_params ] = build_nifti_params(imginfo,nip,fov1,fov2,dim3,dim4)
% build_nifti_params is called by readdicom7 and readdicom7_allscan to
% extract key acquistion parameters from dicom for retention and use in mhd
% and nifti files.
% 20150629  TLChenevert.  Plucked code out of readdicom7 (identical in
% readdicom7_allscan).
% 20150630 TLC
% 20150702 TLC. Seems functional, clean-out deadwood.

% 20150106 Start ...
% Currently nifti_params.orientation used to create nifti will be saved in one structure, but
% these WILL BE IN ERROR FOR MULTI-ANGLE-STACKS SERIES - fortunately a rare scientific need.
% 20150706 TLC Put in more checks against missing entities:
% "MagneticFieldStrength", "ExamDescription", "SeriesDescription".

nifti_params.seriesnumber = imginfo.SeriesNumber;
%nifti_params.seriesdescription = imginfo.SeriesDescription;
nifti_params.seriesdescription = tryGetField(imginfo,'SeriesDescription', 'UNK');

nifti_params.locfirstslc = nip.loc(:,1);
nifti_params.loclastslc  = nip.loc(:,end);
nifti_params.slicectr2ctr = sqrt( sum( (nifti_params.locfirstslc - nifti_params.loclastslc).^2) ) / max([1 dim3-1]); % SpacingBetweenSlices
nifti_params.timefirstphase = nip.td(1,1); % Should be (close to) "td" of first phase, 1st slice in imastor structure.
nifti_params.timelastphase  = nip.td(1,end); % Should be (close to) "td" of last phase, 1st slice in imastor structure.
nifti_params.temporalsampling = (nifti_params.timelastphase - nifti_params.timefirstphase) / max([1 dim4-1]); % Meaningless for non-dynamic series, of course.
nifti_params.orientation  = squeeze(nip.orient(:,1)); % Old nifti & mhd not designed for multi-angle stacks.  Keep this info in nifti_params.allloc & .allorientation.
nifti_params.slicethickness = tryGetField(imginfo,'SliceThickness'); % 20150204 returns "[]" if no field present
nifti_params.mracqtype = tryGetField(imginfo,'MRAcquisitionType'); % 20150204 returns "[]" if no field present
nifti_params.fov = [fov1 fov2];
nifti_params.pixelspacing = tryGetField(imginfo,'PixelSpacing'); % 20150204 returns "[]" if no field present
nifti_params.inplanephencdir = tryGetField(imginfo,'InPlanePhaseEncodingDirection'); % 20150204 returns "[]" if no field present
nifti_params.imgtypstring = tryGetField(imginfo, 'ImageType', 'UNK'); % 20150308
nifti_params.rows = tryGetField(imginfo, 'Rows', 0); % 20150308
nifti_params.columns = tryGetField(imginfo, 'Columns', 0); % 20150308

nifti_params.defname = make3TRID(imginfo,'di'); % 'di' will use de-identification xform.  Only alternative is 'id' that returns date-based id.

% % Need to deal with possible replicate locs in 4D.  Select unique slice locs, but do NOT resort via 'unique':
nifti_params.allloc = nip.loc;
nifti_params.allorientation = nip.orient;

% Institutional, system, and exam demographics:
% nifti_params.institution = imginfo.InstitutionName;
nifti_params.institution = tryGetField(imginfo,'InstitutionName','UNK');
%nifti_params.manufacturer = imginfo.Manufacturer;
nifti_params.manufacturer = tryGetField(imginfo,'Manufacturer','UNK');
%nifti_params.model = imginfo.ManufacturerModelName;
nifti_params.model = tryGetField(imginfo,'ManufacturerModelName','UNK');
%nifti_params.fieldstrength = imginfo.MagneticFieldStrength;
nifti_params.fieldstrength = tryGetField(imginfo,'MagneticFieldStrength',0);
nifti_params.stationname = tryGetField(imginfo,'StationName','UNK');
%nifti_params.serialno = imginfo.DeviceSerialNumber;
nifti_params.serialno = tryGetField(imginfo,'DeviceSerialNumber','UNK');
%nifti_params.swversion = imginfo.SoftwareVersion;
nifti_params.swversion = tryGetField(imginfo,'SoftwareVersion','UNK');
%nifti_params.studydesc = imginfo.StudyDescription;
nifti_params.studydesc = tryGetField(imginfo,'StudyDescription','UNK');

% Possible fourth-dimension arrays for retention as an "FYI" in (mhd/nifti) header comment fields.
nifti_params.td = nip.td;   % [nslc nphase] All slice & time-delay info.
nifti_params.ppd = nip.ppd; % [nslc nphase] All slice & alternative time-delay info.
nifti_params.echotime = nip.echotime; % [1 nphase] 1st slice TEs represtative of all slices
nifti_params.naverages = nip.naverages; % [1 nphase] Navs may vary with bvalue, but 1st slice represtative of all. 
nifti_params.bvalue = nip.bvalue; % [1 nphase] 1st slice represtative of all.
nifti_params.diffdir = nip.diffdir; % [4 nphase] 1st slice represtative of all.
nifti_params.imgtyp = nip.imgtyp; % [1 nphase] 1st slice represtative of all.

end

