function [img,label,fov] = readBMP(varargin)

if nargin && ischar(varargin{1})
    [path,fname,ext] = fileparts(varargin{1});
    fname = [fname,ext];
else
    [fname,path] = uigetfile('*.bmp','Load BMP Images');
end
if fname
    % Find files to load:
    % - assume file name is in form: namestr_~#gstr~_rec####.bmp
    fname = fname(1:(end-8));
    us = strfind(fname,'_');
    k = strfind(fname,'~#');
    namestr = fname(1:us(end-1)-1);
    if ~isempty(k) % Retrospective gating - multiple image sets
        % Select which images to load:
        nbins = str2double(fname(k+4:k+5));
        gstr = cell(1,nbins);
        for i = 1:nbins
            gstr{i} = sprintf('~#%02u%02u~',i-1,nbins);
        end
    else
        gstr = {''};
    end
    % Read relevant parameters from rec.log file:
    fid = fopen(fullfile(path,[namestr,'_',gstr{1},'_','rec.log']));
    if fid>0
        d1 = []; d2 = []; d3 = []; voxsz = []; atend = false;
        while (isempty(d1) || isempty(d2) || isempty(d3) || isempty(voxsz)) && ~atend
            line = fgetl(fid);
            if ischar(line)
                line = strsplit(line,'=');
                switch line{1}
                    case 'Result Image Height (pixels)'
                        d1 = str2double(line{2});
                    case 'Result Image Width (pixels)'
                        d2 = str2double(line{2});
                    case 'Pixel Size (um)'
                        voxsz = str2double(line{2});
                    case 'Sections Count'
                        d3 = str2double(line{2});
                end
            else
                atend = true;
            end
        end
        fclose(fid);

        % Initialize image matrix
        d4 = length(gstr);
        img = zeros(d1,d2,d3,d4);
        fov = [d1,d2,d3] * voxsz/1000;
        label = repmat({'CT'},1,d4);

        % Read in each slice
        hw = waitbar(0, 'Loading Slices ...');
        for ig = 1:d4
            fnames = dir(fullfile(path,[namestr,'_',gstr{ig},'_','rec*.bmp']));
            fnames = {fnames(:).name};
            ind = cellfun(@(x) isnan(str2double(x((end-7):(end-4)))) , fnames);
            fnames(ind) = [];
            n = length(fnames);

            % Read all images into matrix
            for i=1:n
                img(:,:,i,ig) = imread(fullfile(path,fnames{i}));
                waitbar(i/n, hw)
            end
            img(isnan(img)) = 0;
        end
        close(hw);
    end
end