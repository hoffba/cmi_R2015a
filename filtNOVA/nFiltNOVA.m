function img = nFiltNOVA(img,varargin)
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
if nargin
    if nargin==1 % Check in GUI is needed to input parameters
        prompt = {'Radius of moving window cube','Filter passes',...
                  'Initial noise estimate','Trapezoidal plateau (cannot be larger than max extent)',...
                  'Trapezoidal max extent'};
        title = 'NOVA filter parameters';
        default = {'2', '2', '200', '0', '3'};
        answer = str2double(inputdlg(prompt,title,1,default));
        if ~isempty(answer) && ~any(isnan(answer))
            n = answer(1);
            run = answer(2);
            sDev0 = answer(3);
            p = answer(4);
            d = answer(5);
        else
            go = false;
        end
    elseif nargin==6
        n    = varargin{1};
        run  = varargin{2};
        sDev0 = varargin{3};
        p    = varargin{4};
        d    = varargin{5};
    else
        go = false;
    end
    if go % Now check that parameters make sense:
        n = round(n);
        run = round(run);
        go = (n>0) && (run>0) && (d>0) && (p<d) && (sDev0>0);
    end
    if go % Perform the NOVA filter
        t = tic;
        
        % Determine max image value
        mmax = max(img(:));
        [dims(1),dims(2),dims(3),nv] = size(img);
        
        % Pad image for filter window radius
        img = padarray(img,[n,n,n,0],'replicate');
        filtIMG = zeros(dims+2*n);
        Wcum = zeros(dims);
        sDev = ones(dims)*sDev0;
        
        % Skips outer voxels:
        yy = (1:dims(1))+n;
        xx = (1:dims(2))+n;
        zz = (1:dims(3))+n;
        ntot = (2*n+1)^3 * run;
        
        for iv = 1:nv
            % Loop over filter window:
            ct = 0;
            h = waitbar(0, 'Filtering image...');
            for m=1:run

                % Perform a NOVA filter:
                for i = -n:n
                    for j = -n:n
                        for k = -n:n
                            W = kernelNOVA(img(i+yy,j+xx,k+zz,iv)-img(yy,xx,zz,iv),sDev,p,d);
                            Wcum = Wcum + W;
                            filtIMG(yy,xx,zz) = filtIMG(yy,xx,zz) + W.*img(i+yy,j+xx,k+zz,iv);

                            % Increment the waitbar
                            ct = ct + 1;
                            waitbar(ct/ntot,h,...
                                ['Filtering: ',num2str(ct),' of ',num2str(ntot)])
                        end
                    end
                end
                % Normalize by cumulative weights:
                filtIMG(yy,xx,zz) = filtIMG(yy,xx,zz) ./ Wcum;
                filtIMG(filtIMG>mmax) = mmax;
                
                % Calculate local noise:
                if m ~= run % Don't run stDev calculation on last run
                    sDev = localStDev(filtIMG-img(:,:,:,iv),n);
                    if length(hsDev)>1
                        edges = linspace(min(sDev(:)),max(sDev(:)),500);
                        k = histc(sDev(:),edges);
                        hpeak = edges(find(k(3:end)==max(k(3:end)),1,'last') + 2);
                        disp(['Peak StdDev = ',num2str(hpeak)])
                        sDev(sDev>hpeak) = hpeak;
                    else
                        disp('Peak StdDev = 0 (no histogram)')
                    end
                end
                % Update img for next run:
                img(:,:,:,iv) = filtIMG;
                Wcum(:) = 0;
                filtIMG(:) = 0;
            end
            sDev(:) = sDev0;
            fsec = round(toc(t));
            fmin = floor(fsec/60); fsec = mod(fsec,60);
            fhrs = floor(fmin/60); fmin = mod(fmin,60);
            disp(['NOVA filter complete: ',sprintf('%02u:%02u:%02u',fhrs,fmin,fsec)]);
        end
        img = img((n+1):(end-n),(n+1):(end-n),(n+1):(end-n),:);
        close(h);
    else
        img = [];
    end
else
    img = [];
end