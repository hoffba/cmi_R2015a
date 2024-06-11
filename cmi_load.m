function [img,label,fov,orient,info,fnameOut] = cmi_load(imgflag,d,fullname)
% loads images for umich2 program
% Inputs: imflag -   type of load request: 0=mask, 1=image, 2=image appended
%         d -        optional input for dimensions if loading a mask
%         fullname - cell array of file names to load
% Outputs: img - double array of image pixel values
%          label - cell array of strings for each vector
%          fov - spatial extents of the image, [1x3] for 3D images

% Accepted file types:
dtypes = {...                                           Image       Mask
            {'.hdr'},                   'ANALYZE',      true,       true;...   % 1
            {'.mhd'},                   'MHD',          true,       true;...   % 2
            {'.fld'},                   'FLD',          true,       true;...   % 3
            {'.fdf'},                   'FDF',          true,       false;...  % 4
            {'.vff'},                   'VFF',          true,       false;...  % 5
            {'fid','2dseq'},            'BrukerMRI',    true,       true;...   % 6
            {'.log'},                   'BrukerCT',     true,       true;...   % 7
            {'.dcm','DICOMDIR'},        'DICOM',        true,       false;...  % 8
            {'.nii.gz','.nii'},         'NIFTI',        true,       true;...   % 9
            {'.sur','.mrd'},            'MRS',          true,       true;...   % 10
            {'.mask'},                  'MASK',         false,      false;...  % 11
            {'.tif'},                   'TIFF',         true,       true;...   % 12
            {'.jpg'},                   'JPG',          true,       false;...  % 13
            {'.mat'},                   'MAT',          true,       false;...  % 14
            {'.fid'},                   'FID',          true,       false;...  % 15
            {'.vox'},                   'VOX',          true,       false;...  % 16
         };

img = [];
label = {};
fov = [];
orient = [];
info = [];
fnameOut = [];
                
if (nargin == 0)
    imgflag = 1;
end
if (nargin < 2)
    d = [];
end
if (nargin < 3)
    fullname = {};
elseif ~isempty(fullname)
    if ischar(fullname)
        fullname = {fullname};
    end
    ind = ~cellfun(@(x)exist(x,'file'),fullname);
    if any(ind)
        error(['Files not found:\n',repmat('   %s\n',1,nnz(ind))],fullname{ind});
    end
end


% Ask for input
gopt = isempty(fullname);
if gopt
    % Determine filter for UI
    if imgflag % Select an image file
        filter = dtypes([dtypes{:,3}],1:2);
    else % Select a VOI file
        filter = dtypes([dtypes{:,4}],1:2);
    end
    % Determine multi-select option
    if imgflag
        str = 'Image';
        ms_flag = 'on';
    else
        str = 'Mask';
        ms_flag = 'off';
    end
    % User selects file(s)
    filter = [{[filter{:,1}],'All Images'};{{''},'All Files'};filter];
    filter(:,1) = cellfun(@(x)sprintf('*%s;',x{:}),filter(:,1),'UniformOutput',false);
    [fullname,fpath,findex] = uigetfile(filter,['Load ',str,':'],'MultiSelect',ms_flag);
    if fpath
        if ischar(fullname) % single selection
            fullname = {[fpath,fullname]};
        else % multi-selection
            if any((findex-1)==[4,6]) % Saved by slice, don't accept multiple selection
                fullname = fullname(1);
            end
            fullname = cellfun(@(x)fullfile(fpath,x),fullname,'UniformOutput',false);
        end
    else
        fullname = {};
    end
end

if ~isempty(fullname)
    fnameOut = fullname{1};
    if strcmp(fnameOut(end-2:end),'.gz')
        fnameOut(end-2:end) = [];
    end
    nf = length(fullname);
    alabel = cell(1,nf);
    ok = 1;
    for i = 1:nf
        if ok
            % Evaluate format-relevant load function
            [fpath,tname,ext] = fileparts(fullname{i});
            % If GNU zipped, unzip to same folder:
            
            % 20200220: CJG Not sure I want this
%             if strcmp(ext,'.gz')
%                 if ~exist(fullname{i}(1:end-3),'file')
%                     disp('Unzipping ...');
%                     gunzip(fullname{i});
%                 end
%                 fullname{i} = fullname{i}(1:end-3);
%                 [~,~,ext] = fileparts(fullname{i});
%             end

            if any(strcmp(tname,{'2dseq','fid'}))
                ext = tname;
            end
            if isempty(ext) || isfolder(fullname{i})
                datatype = 'DICOM';
            else
                idt = find(cellfun(@(x)contains(fullname{i},x),dtypes(:,1)),1);
                if isempty(idt) && isdicom(fullname{i})
                    idt = 8;
                end
                datatype = dtypes{idt,2};
            end
            if i==1
                if any(strcmp(datatype,{'DICOM','FDF','Bruker'}))
                    fnameOut = fpath;
                end
            end
            [timg,tlabel,tfov,torient,tinfo] = feval(['read' datatype],fullname{i},d);
            [dt(1),dt(2),dt(3),dt(4)] = size(timg);
            % Choose images to load from the 4D set
            if gopt && (dt(4) > 1)
                if imgflag % Loading images
                    [sel,ok] = listdlg('PromptString','Choose single image to load:',...
                                        'ListString',tlabel,'SelectionMode','multiple');
                else % Loading masks
                    if all(dt(1:3)==d(1:3))
                        [sel,ok] = listdlg('PromptString','Choose VOI to load:',...
                                            'ListString',tlabel,'SelectionMode','single');
                    else
                        ok = 0;
                        errordlg('Mask is wrong size for this image!');
                    end
                end
                timg = timg(:,:,:,sel);
                tlabel = tlabel(sel);
            end
            if ok && dt(1)>0
                alabel{i} = tlabel;
                if i==1 % first image loaded initializes the image data
                    img = timg;
                    fov = tfov;
                    orient = torient;
                    info = tinfo;
                    d0 = dt;
                else
                    [d1(1),d1(2),d1(3),d1(4)] = size(timg);
                    if all(d0(1:3) == d1(1:3))
                        img = cat(4,img,timg);
                    else
                        errordlg('Images must be the same 3D size!')
                    end
                end
            end
        end
    end
    label = [alabel{:}];
    % if ~isempty(fov)
    %     fov = fov([2,1,3]); % change to matlab coordinates: [y,x,z]
    % end
end
if (imgflag == 0)
    img = logical(img);
end







