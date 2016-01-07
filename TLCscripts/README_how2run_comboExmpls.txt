To run p-combos

(1) Copy distributed "combo" p-code into your home "Matlab" directory
(or any other directory in your Matlab path).

(2) start Matlab and run the command as if executing regular m-function, e.g.,:

	(a) to split "enhanced" (multi-frame) DICOM, browse inside exam directory at the
            level of "DICOMDIR", and execute:
	>> splitenhd_combo; % then choose destination directory and patient
	    NOTE: this will also automatically generate series catalogue "scaninfo.txt"

	(b) to re-catalogue "legacy" DICOM series, browse to one level above "exam" folder
	    containing separate series folders, and remove all "txt" and "mat" files from
	    series/exam, and execute:
	>> scainfo_combo; % choose exam folder to catalogue
	    NOTE: the resulting "scaninfo.txt" can be viewed to select series of interest
		  for subsequent processing

	(c) to extract specific series imageds and meta-data, execute "readdicom7_combo" 
	   (and choose series):
	>> [imdat,dim1,dim2,dim3,dim4,fov1,fov2,fov3,fov4,tr,unique_te,ppd,flip, dinf, nifti_params] = readdicom7_combo(0,1, 'fp');
	   NOTE: you can then query the output structures and save relevant info using Matlab
 		eg, for DTI series:
	>> idat = getsafield(imdat, 'idata'); % 4D image data
	>> bvals = getsafield(imdat, 'bvalue'); % corresponding b-value array
	>> ddirs = getsafield(imdat, 'diffdir'); % diffusion direction cosines

	(d) to extract and save DTI series info, you could use "dtiXtract" combo:
	>> idat bvals ddirs] = dtiXtract('y'); % type "help" function name for info
