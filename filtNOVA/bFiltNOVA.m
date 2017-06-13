function [img,iname,voxsz] = bFiltNOVA(img,varargin)
% filtIMG = nFiltNOVA(img) : User inputs filter parameters with GUI
%       or
% filtIMG = nFiltNOVA(img,n,run,sDev,p,d)
%
% Performs 3D NOVA filtering for noise reduction
% (Schilham et al., 2006)
% Inputs: n = 3D window radius (e.g. [5x5] : n=2)
%         run = number of NOVA filter passes (default = 2)
%         sDev = initial noise guess
%         p = trapezoidal plateau
%         d = trapezoidal max
% Output: filtIMG = filtered image
 
go = true;
if nargin>0
    if nargin<2
        n = 2;
    else
        n = round(varargin{1});
    end
    D = 2*n+1;
    if nargin<3
        run = 2;
    else
        run = round(varargin{2});
    end
    if nargin<4
        sDev0 = 200;
    else
        sDev0 = varargin{3};
    end
    if nargin<5
        p = 0;
    else
        p = varargin{4};
    end
    if nargin<6
        d = 3;
    else
        d = varargin{5};
    end
    if nargin<7
        iname = 'unidentified';
    else
        iname = varargin{6};
    end
    if nargin<8
        voxsz = ones(1,3);
    else
        voxsz = varargin{7};
    end
    go = (n>0) && (run>0) && (d>0) && (p<d) && (sDev0>0);
    if go % Perform the NOVA filter
        t = tic;
       
        % Determine max image value
        mmax = max(img(:));
        [dims(1),dims(2),dims(3),nv] = size(img);
       
        % Pad image for filter window radius
        img = padarray(img,[n,n,n,0],'replicate');
        filtIMG = zeros(dims+2*n);
        Wcum = zeros(dims);
        sDev = ones(dims);
       
        % Skips outer voxels:
        yy = (1:dims(1))+n;
        xx = (1:dims(2))+n;
        zz = (1:dims(3))+n;
       
        for iv = 1:nv
            sDev(:) = sDev0;
            % Loop over filter window:
            for m=1:run
 
                % Perform a NOVA filter:
                ct = 1;
                for i = -n:n
                    for j = -n:n
                        for k = -n:n
                            if mod(ct,30)==0
                                disp([sprintf('%5.1f',ct/D^3*100),'%']);
                            end
                            W = kernelNOVA(img(i+yy,j+xx,k+zz,iv)-img(yy,xx,zz,iv),sDev,p,d);
                            Wcum = Wcum + W;
                            filtIMG(yy,xx,zz) = filtIMG(yy,xx,zz) + W.*img(i+yy,j+xx,k+zz,iv);
                            ct = ct+1;
                        end
                    end
                end
                % Normalize by cumulative weights:
                filtIMG(yy,xx,zz) = filtIMG(yy,xx,zz) ./ Wcum;
                filtIMG(filtIMG>mmax) = mmax;
 
                % Calculate local noise:
                if m ~= run % Don't run stDev calculation on last run
                    sDev = localStDev(filtIMG-img(:,:,:,iv),n);
                    sDev(sDev==0) = 1;
                    edges = linspace(min(sDev(:)),max(sDev(:)),500);
                    k = histc(sDev(:),edges);
                    hpeak = edges(find(k(3:end)==max(k(3:end)),1,'last') + 2);
                    disp(['Peak StdDev = ',num2str(hpeak)])
                    sDev(sDev>hpeak) = hpeak;
                end
                % Update img for next run:
                img(:,:,:,iv) = filtIMG;
                Wcum(:) = 0;
                filtIMG(:) = 0;
            end
            fsec = round(toc(t));
            fmin = floor(fsec/60); fsec = mod(fsec,60);
            fhrs = floor(fmin/60); fmin = mod(fmin,60);
            disp(['NOVA filter complete: ',sprintf('%02u:%02u:%02u',fhrs,fmin,fsec)]);
        end
        img = img((n+1):(end-n),(n+1):(end-n),(n+1):(end-n),:);
    else
        img = [];
    end
else
    img = [];
end