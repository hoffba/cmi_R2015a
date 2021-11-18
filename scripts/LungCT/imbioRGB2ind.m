function img = imbioRGB2ind(fpath,cmap)

if nargin==1
    cmap = [];
end

% Find DICOM files in folder
fname = dir(fullfile(fpath,'*.dcm'));
fname = {fname.name};

% Load images
for i = 1:numel(fname)
    info(i) = dicominfo(fullfile(fpath,fname{i}));
    img(:,:,:,i) = dicomread(info(i));
end
[d(1),d(2),d(3),d(4)] = size(img);

% Permute to 2D for processing
img = reshape(permute(img,[1,2,4,3]),d(1),[],3);

% Remove grayscale
timg = img - repmat(min(img,[],3),1,1,3);

% Find unique color values
tcmap = double(unique(reshape(timg,[],3),'rows'))/2^info(1).BitDepth;

% Convert to indexed image
timg = rgb2ind(timg,tcmap);

% Map to input cmap
if nargin==2 && size(cmap,2)==3
    img = zeros(size(timg));
    rho = corr(cmap',tcmap');
    for i = 1:size(cmap,1)
        ind = find(rho(i,:)>0.95)-1;
        if numel(ind)==1
            img(timg==ind) = i;
        elseif numel(ind) > 1
            warning('Multiple colormap matches found. Adjust tolerance.');
        end
    end
else
    img = timg;
end

% Reshape to image dimensions
d = [d(1),d(2),d(4)];
img = reshape(img,d);

% Save to file
thk = info(2).ImagePositionPatient(3) - info(1).ImagePositionPatient(3);
fov = [info(1).PixelSpacing',thk] .* d;
orient = info(1).ImageOrientationPatient;
orient = [orient(1:3),orient(4:6),cross(orient(1:3),orient(4:6)),info(1).ImagePositionPatient;0 0 0 1];
fname = fullfile(fileparts(fpath),sprintf('%s_%s.nii.gz',info(1).PatientID,info(1).StudyDate));
saveNIFTI(fname,img,{'ImBio_PRM'},fov,orient);