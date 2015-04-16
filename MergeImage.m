function nImg = MergeImage(varargin)
%Merges the two registered staggered images from odd and even into one
%dicom again, takes either CMIClass inputs or takes two filename inputs
if nargin == 2
    Odd = varargin{1};
    Even  = varargin{2};
    if isa(Odd, 'CMIclass') && isa(Even, 'CMIclass')
    
    else
        [imgOdd, label1, fov1] = readMHD(Odd);
        [imgEven, ~, ~] = readMHD(Even);
        dims = size(imgOdd);
        nImg = zeros(dims(1), dims(2), 2*dims(3));
        %for loop route, reshape may be faster but as there are few slices
        %it doesn't matter in this case
        for x = 1:dims(3)
            c = 2*x;
            nImg(:,:,c-1) = imgOdd(:,:,x);
            nImg(:,:,c) = imgEven(:,:,x);
        end
        [fname, pname] = uiputfile('*.mhd', 'Save New File Where?', label1{1});
        saveMHD(fullfile(pname, fname), nImg, dims, fov1);
    
    end
elseif nargin == 1
    cmi = varargin{1};
    img = cmi.img.mat;
    imgOdd = img(:,:,:,1);
    imgEven = img(:,:,:,2);
    dims = size(imgOdd);
    nImg = zeros(dims(1), dims(2), 2*dims(3));
    fov1 = cmi.img.voxsz .* dims;
    label1 = cmi.img.name;
    for x = 1:dims(3)
        c = 2*x;
        nImg(:,:,c-1) = imgOdd(:,:,x);
        nImg(:,:,c) = imgEven(:,:,x);
    end
    [fname, pname] = uiputfile('*.mhd', 'Save New File Where?', label1);
    saveMHD(fullfile(pname, fname), nImg, dims, fov1);
    
    
end

end

