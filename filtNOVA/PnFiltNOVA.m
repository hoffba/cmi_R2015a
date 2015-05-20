function img = PnFiltNOVA(img,varargin)
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
        prompt = {'Radius of moving window cube',...
                  'Filter passes','Initial noise estimate',...
                  'Trapezoidal plateau (cannot be larger than max extent)',...
                  'Trapezoidal max extent'};
        title = 'NOVA filter parameters';
        default = {num2str(1:size(img,4)),'2', '2', '200', '0', '3'};
        answer = inputdlg(prompt,title,1,default);
        if ~isempty(answer) && ~any(isnan(str2double(answer(2:end))))
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
        hw = waitbar(0,'Prepping data for parallel processing ...');
        [dims(1),dims(2),dims(3),nv] = size(img);
        nstep =  run * nv;
        
        % Prep data for parallel processing
        mmax = max(img(:));
%         nwrk = min(floor(dims(3)/(4*n+2)),12);
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            poolobj = parpool;
        end
        nwrk = poolobj.NumWorkers;
        sDev = distributed.ones(dims(1:3));
        dimg = distributed.ones(dims(1:3));
        for iv = 1:nv
            sDev(:) = sDev0;
            dimg(:) = img(:,:,:,iv);

            % Loop over filter iterations
            for m = 1:run
                waitbar((m+run*(iv-1)-1)/nstep,hw,'Performing filter on distributed matrix ...');
                disp(['Run #',num2str(m)]);
                spmd(nwrk)
                    % Initialize necessary variables
                    labp = mod(labindex,numlabs)+1;
                    labm = mod(labindex-2,numlabs)+1;
                    locimg = getLocalPart(dimg);
                    locsDev = getLocalPart(sDev);
                    mdata = labSendReceive(labp,labm,locimg(:,:,(end-n+1):end));
                    pdata = labSendReceive(labm,labp,locimg(:,:,1:n));
                    if labindex==1
                        mdata = repmat(locimg(:,:,1),[1,1,n]);
                    elseif labindex==numlabs
                        pdata = repmat(locimg(:,:,end),[1,1,n]);
                    end
                    padimg = padarray(cat(3,mdata,locimg,pdata),[n,n,0],'replicate');
                    locdims = size(padimg);
                    filtImg = zeros(locdims);%padimg;
                    yy = (n+1:locdims(1)-n);
                    xx = (n+1:locdims(2)-n);
                    zz = (n+1:locdims(3)-n);

                    % Perform a NOVA filter:
                    Wcum = zeros(locdims-2*n);
                    for i = -n:n
                        for j = -n:n
                            for k = -n:n
                                W = kernelNOVA(padimg(i+yy,j+xx,k+zz)-padimg(yy,xx,zz),...
                                               locsDev,p,d);
                                Wcum = Wcum + W;
                                filtImg(yy,xx,zz) = filtImg(yy,xx,zz) ...
                                                    + W.*padimg(i+yy,j+xx,k+zz);
                            end
                        end
                    end

                    % Normalize by cumulative weights:
                    filtImg(yy,xx,zz) = filtImg(yy,xx,zz) ./ Wcum;
                    filtImg(filtImg>mmax) = mmax;

                    % Update the codistributed array
                    codist = getCodistributor(dimg);
                    a = sum(codist.Partition(1:labindex-1)) + 1;
                    b = a + codist.Partition(labindex) - 1;
                    dimg(:,:,a:b) = filtImg(n+1:end-n,n+1:end-n,n+1:end-n);
                    
                    % Calculate local noise:
                    if m ~= run % Don't run stDev calculation on last run
                        locsDev = localStDev(filtImg-padimg,n);
                        edges = linspace(min(locsDev(:)),max(locsDev(:)),500);
                        k = histc(locsDev(:),edges);
                        hpeak = edges(find(k(3:end)==max(k(3:end)),1,'last') + 2);
                        disp(['Peak StdDev = ',num2str(hpeak)])
                        locsDev(locsDev>hpeak) = hpeak;
                        sDev(:,:,a:b) = locsDev;
                    end
                end
            end

            img(:,:,:,iv) = gather(dimg);

            disp(['NOVA filter complete: ',datestr(toc(t)/(60*60*24),'HH:MM:SS')]);
            
%             fsec = round(toc(t));
%             fmin = floor(fsec/60); fsec = mod(fsec,60);
%             fhrs = floor(fmin/60); fmin = mod(fmin,60);
%             disp(['NOVA filter complete: ',sprintf('%02u:%02u:%02u',fhrs,fmin,fsec)]);
        end
        delete(hw);
    else
        img = [];
    end
else
    img = [];
end