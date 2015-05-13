function stat = batch_hetProc(img,mask,fov,fnout,varargin)
% Add batch job to process Minkowski Functionals and Entropy
% Inputs:   img   = 3D image to process
%           mask  = binary image mask
%           fov   = image extents
%           fnout = output image file name

[fpath,fname] = fileparts(fnout);
if nargin<5
    % This code starts the batch
    disp(['Starting batch heterogeneity analysis process ... ',fname])
    batch(@batch_hetProc,1,{img,mask,fov,fnout,1});
    disp(' ... done')
elseif varargin{1}==1
    % This code is run inside the batch
    disp(fname);
    R = 4;

    disp('Calculating Minkowski Functionals ...')
    nMF = (R*2+1)*ones(1,3);
    thresh = -950:100:-300;
    [fimg,labels] = minkowskiFun(img,nMF,thresh,mask);
    
    d = size(fimg);
    stat = false(d(4:5));
    
    disp('Saving Minkowski Functional Maps ...')
    for i = 1:size(fimg,5)
        for j = 1:size(fimg,4)
            oname = [fname,'_t',num2str(thresh(j)),'_',labels{i},'.mhd'];
            stat(j,i) = saveMHD(fullfile(fpath,oname),fimg(:,:,:,j,i),'',fov);
        end
    end

    disp('Calculating Local Entropy ...')
    nE = bwellipsoid(R*ones(1,3));
    fimg = entropyfilt(img,nE);
    scl = 32000/max(fimg(:));
    fimg = fimg*scl;
    disp(['Entropy scaled by: ',num2str(scl)])

    disp('Saving Entropy Results ...')
    saveMHD(fullfile(fpath,[fname,'_Entropy.mhd']),fimg,[],fov);
end
