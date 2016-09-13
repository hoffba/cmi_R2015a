function examoutdir = split_series(examindir)
%
% TLC 20140807: Spun-off of show_scan_info.  Use to convert single folder
%               containing all images from several series into distinct series folders.
%               Add flag to auto run "show_scan_info"
%               Swith to run "read_enh" if the first IM* file found to be ENHANCED DICOM.
%
% MD 20140820: "qiba_split" called by "qiba_show_scan" only for "simple" dicom
% MD 20140826: replaced the code-block below with a function
% MD 20140826: allowed independent script execution, but pay attention to
% hard-coded options below
% MD & TLC 20140902:    Taken from QIBA SW Project
% TLC 20141223:  Use wildcard = 'ss' as relay to dcmfsearch.m to NOT apply
% filesize-based filter to reject inconsistent dicom entities from given series folder.

global fname dtes examID lastdcdir; % This is so I can label images with appropriate filenames.

startdir = pwd;
% CHANGE wildwildcard to 'y' IN BOTH THIS SCRIPT AND readdicom7_allscan.
wildwildcard = 'n'; % 20140502. def = 'n'.  'y' means image filenames found by wide-open "*" wildcard.  Be sure to manually delete all extraneous files (infos, jpgs, mats,...).
% wildwildcard = 'ss'; % used to be 'n'.  FOR NOW, NOT DOING THIS 'ss' fliter 20141223 
run_show_scan_info = 'n'; % 'y' = also run show_scan_info.
run_enh = 'y'; % MD 20140826: will not run "enhanced" for "qiba"

% if (nargin < 1) % MD 20140822: run from "qiba_show_scan_info"
%     examindir = startdir;
% end


%select exam of interest
if (nargin == 0)       % No input argument, MD 20140826: run by itself
    browz = 0; % Assign browse starting point to default if browz is empty
    disp('  ***   Select EXAM FOLDER containing all series folders and/or images ***');
    disp('  (FYI) Do NOT go so deep as to select a SERIES FOLDER! ');
    examindir = uigetdir(pwd, 'Pick Exam of Interest');
    %examindir = uigetdir('C:\data', 'Pick Exam of Interest');
else
%     examindir = lastdcdir;
    examindir = startdir; % MD 20140826: run from "qiba_show_scan_info"
end

cd(examindir);
[basepath, basefolder, baseext] = fileparts(pwd);
cd ..
%[basepath, basefolder, baseext] = fileparts(pwd);
basedir = pwd;

% Define destination exam folder name
examfolder = [basefolder '_SeriesSplit'];

% Check if examfolder already exists, if not make it.
if (isdir(examfolder))
    disp(['  (FYI) Using existing split-series folder = "' examfolder '"'])
else
    mkdir(examfolder);
end % if isdir

cd(examfolder);
examoutdir = pwd;
%[examoutpath,examoutfolder,examoutext] = fileparts(pwd);
scanlist = dir();
%scanlist = scanlist(3:end, :);   % remove the directories "." and ".." (assumes "Windows")
chkdir = getsafield(scanlist, 'isdir')'; % MD 20140820
% MD 20140820: check if split series exist
if (length(chkdir(chkdir == 1)) > 3) % split series exist in the folder
    disp(['  (FYI) Found split series to use in = "' examoutdir '"'])
else

% Go back into examindir and start parse process
    cd(examindir);
imagelist = dcmfsearch(wildwildcard); % MD 20140826: replaced the code-block below with a function
%imagelist = dir('I*0*'); % PENDING - Start with limited template that works for now.  Generalize later.
%     imagelist = dir('*.dcm');
%     if ( isempty(imagelist) )
%     % Choose from next two lines that detects dicom images w/o counting image_order.mat file.
%         imagelist = dir('I*0*');  % Paris has I000## images.  This will catch both HFH and Paris.
%         if ( isempty(imagelist) )
%         % 20140502 Workaround when images fnames start with multiple integers (like 001, ..., 100, 101,...).
%         % Then use a simple " * " wildcard:
%             if ( wildwildcard == 'y' )
%                 junk = dir();
%                 imagelist = junk(3:end, :); % Need to skip over "." and ".." included in list.
%             else       
%                 for jjj = 0:9 % Then count files that begin with an integer.  Once found, go no further
%                     if ( isempty(imagelist) ) 
%                         wildcard = [num2str(jjj) '*'];
%                         imagelist = dir(wildcard);
%                     end % if isempty
%                 end % for jjj
%             end % if wildwildcard
%         end % if isempty
%     end % if isempty
        if (isempty(imagelist)) % above wild cards did not work
            disp(['ERROR: empty DICOM image file list for "' dirname '" directory ... QUITTING!']) 
            fclose(fid);% error out gracefully
            return;
        end
    numimgs = length(imagelist);
    if (numimgs == 0)
        disp(' ERROR: No image data in selected exam folder. QUITTING! ');
        return;
    end

% Let's check to see if these images are Enhanced DICOM
    di1 = dicominfo(imagelist(1).name);

    if (isfield(di1,'PerFrameFunctionalGroupsSequence'))
        % This appears to be Enhanced.  Jump up 1 directory level and run "read_enh(4)".
        if (run_enh == 'y') % allow "enhanced" input
            cd ..;
            read_enh(4);     % Call read_enh using newly-created "diropt = 4"
            serieslist = []; % Not used, but assign anway
            folderindex = 0; % Not used, but assign anway
            noserno = 0;     % Not used, but assign anway
        else % "enhanced" is not allowed
            disp(' ERROR: image data in compressed DICOM format. QUITTING!');
            disp('  (FYI) uncompressed DICOM required for processing'); 
            return;
        end
   
    else % Then probably all images are non-enhanced DICOM - assume so...
    
        folderindex = 0; % Count folders 1,2,3...
    % %loop thru each scan and display scan info
        serieslist = [];
        noserno = 0;
        for ii = 1:numimgs

            if ( imagelist(ii).isdir == 1) % 20140331
                folderindex = folderindex + 1; % Keep count of directories, but otherwise skip over them.

            else
                di = dicominfo(imagelist(ii).name);
                if (isfield(di,'SeriesNumber'))
                    serno = di.SeriesNumber;
                    serfoldercheck = [examoutdir '\' num2str(serno)];
                    if (isdir(serfoldercheck))
                        copyfile(imagelist(ii).name,serfoldercheck);
                    else
                        mkdir(examoutdir,num2str(serno));
                        serieslist = [serieslist; serno]; % Add new series to list
                        copyfile(imagelist(ii).name,serfoldercheck);
                    end % if isdir          
                else
                    noserno = noserno + 1; % Keep count of files without series number.
                end % if isfield
            end % if imagelist(ii).isdir == 1
        end % for ii
%     disp(['  (FYI) number of dicom images that do not have a SeriesNumber field = "' num2str(noserno) '"']);
%     disp(['  (FYI) number of sub-folders skipped over = "' num2str(folderindex) '"']);
    %disp(' ');
        disp(['  (FYI) Split exam series saved in = "' examoutdir '"'])
    %disp(' ');

    if (run_show_scan_info == 'y')
        lastdcdir = examoutdir;
        show_scan_info(1);
    end % if run_show_scan_info
    
    end % if isfield(di1,'PerFrameFunctionalGroupsSequence')

end % split series exist

cd(startdir);

return;