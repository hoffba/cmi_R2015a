function [img,label,fov] = readBruker(varargin)

go = false;
img = [];
label = {};
fov = [];

% First determine basename:
fname = varargin{1};
[path,fname,ext] = fileparts(fname);
if ~isnan(str2double(fname(end-3:end)))
    go = true;
    recchk = ~isempty(strfind(fname,'_rec'));
    bname = fname(1:end-4-recchk*4);
    serchk = false; nser = 1;
    if recchk
        ind = strfind(fname,'~#');
        serchk = ~isempty(ind);
        if serchk
            nser = str2double(fname((1:2)+ind+3));
            bname = bname(1:end-7);
        end
    end
end

% Next find relevant files in the folder:
if go
    searchstr = bname;
    if serchk
        iserstr = length(bname)+(3:4);
        searchstr = [searchstr,'~#00',sprintf('%02u',nser),'~_rec'];
        nserstr = mat2cell(num2str((1:nser)','%02u'),ones(1,nser));
    elseif recchk
        searchstr = [searchstr,'_rec'];
    end
    fnames = cell(1,nser);
    keep = false(1,nser);
    for iser = 1:nser
        if serchk
            searchstr(iserstr) = nserstr{iser};
        end
        tnames = dir(fullfile(path,[searchstr,'*',ext]));
        if ~isempty(tnames)
            keep(iser) = true;
            tnames = {tnames(:).name};
            if ~recchk
                rnames = dir(fullfile(path,[searchstr,'*_rec*',ext]));
                tnames = setdiff(tnames,{rnames(:).name});
            end
            ind = cellfun(@(x)isnan(str2double(x(end-7:end-4))),tnames);
            fnames{iser} = tnames(~ind);
        end
    end
    fnames = fnames(keep);
    if serchk
        nserstr(~keep) = [];
    end
end

% Read image parameters from .log file:
if go
    % Set default values:
    dims = [zeros(1,2),length(fnames{1}),length(fnames)];
    
    % Try to read log file:
    strs = {'[Reconstruction]'              ,'[Acquisition]';...
            'Pixel Size (um)'               ,'Image Pixel Size (um)';
            'Result Image Width (pixels)'   ,'Number Of Columns';
            'Result Image Height (pixels)'  ,'Number Of Rows'};
    imtype = 2;
    label = 'Proj';
    if recchk
        imtype = 1;
        label = 'CT';
    end
    if serchk
        label = strcat(label,'#',nserstr);
    else
        label = {label};
    end
    fid = fopen(fullfile(path,[fname(1:end-4),'.log']));
    if fid>0
        rchk = false;
        done = false;
        str = fgetl(fid);
        while ischar(str)
            if ~rchk
                rchk = strcmp(str,strs{1,imtype});
            elseif strcmp(str(1),'[')
                done = true;
            else
                str = strsplit(str,'=');
                if strcmp(str{1},strs{2,imtype})
                    voxsz = str2double(str{2});
                elseif strcmp(str{1},strs{3,imtype})
                    dims(2) = str2double(str{2});
                elseif strcmp(str{1},strs{4,imtype})
                    dims(1) = str2double(str{2});
                end
            end
            str = fgetl(fid);
            done = done || ~ischar(str);
        end
        fclose(fid);
    else
        img = imread(fullfile(path,[fname,ext]));
        dims(1:2) = size(img);
    end
    
    % Load all image files:
    hw = waitbar(0,'');
    ntot = prod(dims(3:4));
    img = zeros(dims);
    for i = 1:length(fnames)
        nf = length(fnames{i});
        if nf==dims(3) % Check that # slices matches
            for islc = 1:nf
                img(:,:,islc,i) = imread(fullfile(path,fnames{i}{islc}));
                waitbar((dims(3)*(i-1)+islc)/ntot,hw,...
                    ['Vec #',num2str(i),' ; Slc #',num2str(islc)]);
            end
        end
    end
    close(hw);
    fov = dims(1:3)*voxsz/1000;
end


