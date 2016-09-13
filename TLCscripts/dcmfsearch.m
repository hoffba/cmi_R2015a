function filelist = dcmfsearch(wildwildcard) 
%
% searches for DCM files by wildcard in the current dir
% MD 20141223. Filter-out "not-an-image" DCM entities 
% TLC 20150218.  Extend grope thru wildcards to include *.IMA from Siemens.
% MD 20150226: "exclude" files by mask logic
% TLC 20150611: Invisible Thumbs.db causing problems.  Add to exclusion list.
% MD 20160520: extend "wrong-size" image filer to exclude (by ImageType)
% LOCALIZER/PROJECTION images (inserted in Siemens series)
% TLC 20160629: add *DS_Store' to the exclusion list.

xmsk = {'*bmp' '*jpg' '*png' '*gif' '*tiff' '*txt' '*pdf' '*.db' '*nfo' '*DS_Store'}; % file masks to exclude  

flx = dir('*mat'); % init file list to exclude ".mat" files

for ix = 1:length(xmsk)
    flx = [flx; dir(xmsk{ix})];
end

flx1 = struct2cell(flx); % cast to cell array
% fnmx = {};
% for ixx = 1:length(flx)
%     fnmx{ixx} = flx(ixx).name; % file names to exclude
% end

flall = dir('*'); % all files
fla1 = struct2cell(flall); % cast to cell array

isd = getsafield(flall, 'isdir'); % directories
kd = find(isd>0); % directory entry index

lenfla = length(flall); % length of all-file-list

% fnma = {};
% for i=1:lenfla    
%     fnma{i}=flall(i).name; % all file names
% end
% isx = ismember(fnma,fnmx);
isx = ismember(fla1(1,:),flx1(1,:)); % find the file names to exclude
kx = find(isx > 0); % exclude file index (in the all-file-list)

kin = setdiff(1:lenfla,union(kd,kx));

if (~isempty(kin))
    filelist = flall(kin); % filtered file-list
else
    filelist = [];
end

clear flx flall flx1 fla1 kd kx kin isd lenfla %fnmx fnma
% filelist = dir('*.dcm');
% if ( isempty(filelist) )
%     filelist = dir('IM_*'); % Want to avoid counting "image_order.mat" as a dicom file
%     if ( isempty(filelist) )
%         % filelist = dir('IM0*'); % Want to avoid counting "image_order.mat" as a dicom file
%         filelist = dir('I*0*'); % Want to avoid counting "image_order.mat" as a dicom file
%         if ( isempty(filelist) )
%             
% 
%             filelist = dir('*.IMA'); % 20150218 *IMA. This avoids counting image_order.mat
%             
%             if ( isempty(filelist) ) % 20150218
%             
%                 if (wildwildcard == 'y')
%                     junk = dir();
%                     filelist = junk(3:end, :); % to skip over "." and ".."
%                 else
%                     for jjj = 0:9 % Then count files that begin with an integer.  Once found, go no further
%                         if ( isempty(filelist) ) 
%                             wildcard = [num2str(jjj) '*'];
%                             filelist = dir(wildcard);
%                         end % if isempty
%                     end % for jjj
%                 end % if wildwildcard
%                 
%             end % if isempty % 20150218
%         end % if isempty
%     end % if isempty
% end % if isempty

% % % 20141223 Start ...
if (~isempty(filelist))

fn1 = filelist(1).name;
fnl = filelist(end).name;
inf1 = dicominfo(fn1); infl=dicominfo(fnl);
sameser = ((inf1.SeriesNumber - infl.SeriesNumber)==0);
%sameser
if (sameser == 1) % do check only for single series "folder"
    % disp('single series folder')
    kbtsz = fix(getsafield(filelist,'bytes')/1000); %image size in Kb
    lnfl = length(filelist);
    %kdsz = find(abs(kbtsz-mode(kbtsz))>2*std(kbtsz));
    %%% assuming that "images" are more frequent than non-images
    kdsz = find(abs(kbtsz-mode(kbtsz))>3); % 2^11 - byte-size difference (assuming 12-bit image depth)
    %kdsz  % debug
    knim = [0]; % init "non-image" index
    if ~isempty(kdsz)
        chk1 = dicominfo(filelist(kdsz(1)).name); % dicom for sure
        % MD 20160520: exclude LOCALIZER/PROJECTION images (inserted in Siemens series)
        % if (~isempty(chk1.BitDepth)) % 20160520 "outlier" is the legit "image"
        if (~isempty(chk1.BitDepth) && isempty(strfind(chk1.ImageType,'LOCALI')) && isempty(strfind(chk1.ImageType,'PROJE'))) % 20160520 "outlier" is the legit "image"
            kdsz = setdiff(1:lnfl, kdsz); % flip the index for "non-image" suspects
        end
        %kdsz  % debug
        if (isempty(kdsz)) kdsz = 1; end; % just in case
        for ik = 1:length(kdsz)
            %if (isdicom(filelist(kdsz(ik)).name))
                %chkim = dicomread(filelist(kdsz(ik)).name);
                chkim = dicominfo(filelist(kdsz(ik)).name); % dicom for sure
                %if (prod(size(chkim))<1)
                % MD 20160520: exclude LOCALIZER/PROJECTION images (inserted in Siemens series)
                % if (isempty(chkim.BitDepth)) % 20160520
                if (isempty(chkim.BitDepth) || ~isempty(strfind(chkim.ImageType,'LOCALI')) || ~isempty(strfind(chkim.ImageType,'PROJE'))) % 20160520
                    knim = [knim kdsz(ik)];
                end
            %else
            %    knim = [knim kdsz(ik)];
            %end        
        end    
    else % all images of the same size 
        chk1 = dicominfo(filelist(1).name); % dicom for sure
        if (isempty(chk1.BitDepth)) % if first file is "non-image"
            knim = [knim 1:lnfl]; % all files are non-legit
        end % otherwise all files are legit "images"
    end
    %knim  %debug
    if (length(knim)>1) % remove the non-images from the list
        newk = setdiff(1:lnfl,knim(2:end));
        if (~isempty(newk))
            filelist = filelist(newk); % updated filelist
        else % all files are "non-images"
            filelist = []; % return empty list
            disp('ERROR: No DICOM images found in the folder  ');
            disp(pwd);
        end
    end
end

else % empty folder
    filelist=[];
    disp('ERROR: empty folder  ');
    disp(pwd);
end
% filelist = newlist;
% % % ... End 20141223.

return;