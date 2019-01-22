function [img,label,fov,fnameOut,info] = cmi_load(imgflag,d,fullname)
% loads images for umich2 program
% Inputs: imflag -   type of load request: 0=mask, 1=image, 2=image appended
%         d -        optional input for dimensions if loading a mask
%         fullname - cell array of file names to load
% Outputs: img - double array of image pixel values
%          label - cell array of strings for each vector
%          fov - spatial extents of the image, [1x3] for 3D images
%          fnameOut
%          info - structure containing image metadata


info = [];
% Accepted file types:
if imgflag
    str = 'Image';
else
    str = 'Mask';
end
dtypes = {...
            {'.hdr'},                   'ANALYZE';...   % 1
            {'.mhd'},                   'MHD';...       % 2
            {'.fld'},                   'FLD';...       % 3
            {'.fdf'},                   'FDF';...       % 4
            {'.vff'},                   'VFF';...       % 5
            {'fid','2dseq'},            'BrukerMRI';... % 6
            {'.log'},                   'BrukerCT';...  % 7
            {'.dcm','.1','DICOMDIR'},   'DICOM';...     % 8
            {'.nii'},                   'NIFTI';...     % 9
            {'.sur','.mrd'},            'MRS';...       % 10
            {'.mask'},                  'MASK';...      % 11
            {'.tif'},                   'TIFF';...      % 12
            {'.jpg'},                   'JPG';...       % 13
            {'.mat'},                   'MAT';...       % 14
            {'.fid'},                   'FID';...       % 15
            {'.vox'},                   'VOX';...       % 16
         };
imggp = true(13,1); imggp(10) = false;
maskgp = false(13,1); maskgp([1:3,6,7,9,10,12]) = true;

if (nargin == 0)
    imgflag = 1;
end
if (nargin < 2)
    d = [];
end
if (nargin < 3) || (~iscell(fullname))
    fullname = {};
end

% Check that all input file names exist
gopt = true;
if ~isempty(fullname)
    gopt = any(~cellfun(@(x)exist(x,'file'),fullname));
end

% If not a valid fname, ask for input
if gopt
    % Determine filter for UI
    if imgflag % Select an image file
        filter = dtypes(imggp,:);
    else % Select a VOI file
        filter = dtypes(maskgp,:);
    end
    % Determine multi-select option
    if imgflag
        str = 'on';
    else
        str = 'off';
    end
    % User selects file(s)
    filter = [{[filter{:,1},'.gz'],'All Images'};{{''},'All Files'};filter];
    filter(:,1) = cellfun(@(x)sprintf('*%s;',x{:}),filter(:,1),'UniformOutput',false);
%     [fullname,fpath,findex] = uigetfile([{'*','All'};filter],['Load ',str,':'],'MultiSelect',str);
    [fullname,fpath,findex] = uigetfile(filter,['Load ',str,':'],'MultiSelect',str);
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
            if strcmp(ext,'.gz')
                if ~exist(fullname{i}(1:end-3),'file')
                    disp('Unzipping ...');
                    gunzip(fullname{i});
                end
                fullname{i} = fullname{i}(1:end-3);
                [~,~,ext] = fileparts(fullname{i});
            end
            if any(strcmp(tname,{'2dseq','fid'}))
                ext = tname;
            end
            if isempty(ext)
                datatype = 'DICOM';
            else
                idt = find(cellfun(@(x)any(strcmpi(ext,x)),dtypes(:,1)),1);
                if isempty(idt) && isdicom(fullname{i})
                    idt = 8;
                end
                datatype = dtypes{idt,2};
            end
            if i==1
                if any(strcmp(datatype,{'DICOM','FDF','Bruker'}))
                    fnameOut = fpath;
                end
                cd(fpath);
            end
            [timg,tlabel,tfov,info] = feval(['read' datatype],fullname{i},d);
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
            if ok
                alabel{i} = tlabel;
                if i==1 % first image loaded initializes the image data
                    img = timg;
                    fov = tfov;
                    d0 = dt;
                else
                    [d1(1),d1(2),d1(3),d1(4)] = size(timg);
                    if all(d0(1:3) == d1(1:3))
                        img = cat(4,img,timg);
                    else
                        errordlg('Images must be the same 3D size!')
                    end
                end
            else
                img = [];
                alabel = {};
                fov = [];
                fnameOut = [];
            end
        end
    end
    label = [alabel{:}];
else
    img = [];
    label = {};
    fov = [];
    fnameOut = [];
end
if (imgflag == 0)
    img = img>min(img(:));%logical(img);
end







