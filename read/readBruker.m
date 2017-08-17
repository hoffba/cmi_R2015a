function [img,label,fov] = readBruker(varargin)

go = false;
img = [];
label = {'CT'};
fov = [];
zstep = 1;
dmult = [1,1];

% First determine basename and find all related files:
fname = varargin{1};
[fpath,bname,~] = fileparts(fname);
% fnames = dir(fullfile(path,[bname,'*',ext]));

% Read image parameters from .log file:
fid = fopen(fullfile(fpath,[bname,'.log']));
if fid>2
    info = strsplit(fread(fid,inf,'*char')','\n');
    fclose(fid);
    ind = find(strncmp('[Reconstruction]',info,16),1);
    if isempty(ind) || ~strcmp(bname(end-2:end),'rec') % Acquisition -- for projections
        str2 = strsplit(info{strncmp('Image Pixel Size',info,16)},'=');
        voxsz = str2double(str2{2});
        str2 = strsplit(info{strncmpi('Number Of Rows',info,14)},'=');
        dims(1) = str2double(str2{2});
        str2 = strsplit(info{strncmpi('Number Of Columns',info,17)},'=');
        dims(2) = str2double(str2{2});
        str2 = strsplit(info{strncmpi('Number Of Files',info,15)},'=');
        dims(3) = str2double(str2{2});
        ext = '.tif';
        istart = 0;
    else % Reconstruction -- for reconstructed images
        info(1:ind-1) = [];
        str2 = strsplit(info{strncmp('Pixel Size',info,10)},'=');
        voxsz = str2double(str2{2});
        str2 = strsplit(info{strncmp('Result Image Height',info,19)},'=');
        dims(1) = str2double(str2{2});
        str2 = strsplit(info{strncmp('Result Image Width',info,18)},'=');
        dims(2) = str2double(str2{2});
        str2 = strsplit(info{strncmp('Sections Count',info,14)},'=');
        dims(3) = str2double(str2{2});
        str2 = strsplit(info{strncmp('First Section',info,13)},'=');
        istart = str2double(str2{2});
        str2 = strsplit(info{strncmp('Section to Section Step',info,23)},'=');
        zstep = str2double(str2{2});
        str2 = strsplit(info{strncmp('Filename Index Length',info,21)},'=');
        fnind = strtrim(str2{2});
        str2 = strsplit(info{strncmp('Result File Type',info,16)},'=');
        ext = ['.',lower(strtrim(str2{2}))];
        str2 = strsplit(info{strncmp('Undersampling factor',info,20)},'=');
        usf = str2double(str2{2});
    end
    go = true;
else
%     voxsz = 1;
%     dims = size(imread(fullfile(path,fname)));
end

% Load all image files:
if go
    % Load images by what's in the folder
    fnames = dir(fullfile(fpath,[bname,'*',ext]));
    fnames = {fnames(:).name}';
    fnames(cellfun(@(x)isempty(regexp(x,'_rec\d\d\d\d','once')),fnames)) = [];
    dims(3) = length(fnames);
    
    hw = waitbar(0,'');
    set(findall(hw,'type','text'),'Interpreter','none');
    img = zeros(dims);
    for i = 1:dims(3)
        img(:,:,i) = imread(fullfile(fpath,fnames{i}));
%         img(:,:,i) = imread(fullfile(fpath,sprintf([bname,'%0',fnind,'u',ext],...
%                                                   istart/usf+(i-1)*zstep)));
        waitbar(i/dims(3),hw,{bname;['Loading slice:',num2str(i),' of ',...
                              num2str(dims(3))]});
    end
    close(hw);
    fov = dims.*[dmult,zstep]*voxsz/1000;
    fov = fov([2,1,3]);
end


