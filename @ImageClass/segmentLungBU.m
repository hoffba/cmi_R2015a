% ImageClass function
function segmentLung(self,vec)
pcheck = false; % option to pause at certain points for checking VOI
% Automated lung segmentation algorithm
if any(self.mask.mat(:)) % Must first have a VOI on 2 major brochi near lung
    prompt = {'Custom','Clin. Insp.','Clin. Exp.','Mouse Insp.','Mouse Exp.'};
    [sel,ok] = listdlg('ListString',prompt,'InitialValue',5,'SelectionMode','single');
    if ok
        switch sel
            case 1 % Custom - inputdlg
                sm = 3;
                b = self.scaleB(vec);
                Ti = 607 + b;
                Tt = 125 + b;
                bt1 = 8; % Trachea dilation
                bt2 = 8; % Lung erosion
                bt3 = 9; % End dilation
                bt4 = 3; % End erosion
                slc = 175;%round(self.dims(3)/2);
                prompt = {'Image Smoothing:',...
                          'Lung Threshold:',...
                          'Trachea Threshold:',...
                          'Trachea Smoothing (0-1):',...
                          'Lung Erosion (0-1):',...
                          'End Dilation (0-1):',...
                          'End Erosion (0-1):',...
                          'Display (slice#):'};
                defs = {num2str(sm),...
                        'auto',...%num2str(Ti),...
                        'auto',...%num2str(Tt),...
                        num2str(bt1),...
                        num2str(bt2),...
                        num2str(bt3),...
                        num2str(bt4),...
                        num2str(slc)};
                answer = str2double(inputdlg(prompt,'Lung Segmentation Options',1,defs));
                for i = 1:length(answer) % parse inputs
                    val = answer(i);
                    nanchk = isnan(val);
                    switch i
                        case 1
                            if ~nanchk
                                sm = val;
                            end
                        case 2
                            if nanchk
                                Ti = [];
                            else
                                Ti = val;
                            end
                        case 3
                            if nanchk
                                Tt = [];
                            else
                                Tt = val;
                            end
                        case 4
                            if ~nanchk && (val>0)
                                bt1 = val; % dilation
                            end
                        case 5
                            if ~nanchk && (val>0)
                                bt2 = val; % erosion
                            end
                        case 6
                            if ~nanchk && (val>0)
                                bt3 = val; % dilation
                            end
                        case 7
                            if ~nanchk && (val>0)
                                bt4 = val; % erosion
                            end
                        case 8
                            val = round(val);
                            if ~nanchk && (val>0) && (val<=self.dims(3))
                                slc = val;
                            else
                                slc = [];
                            end
                    end
                end
            case 2 % Clin. Insp.
                sm = 3;
                b = self.scaleB(vec);
                Ti = 607 + b;
                Tt = 125 + b;
                bt1 = 8; % Trachea dilation
                bt2 = 8; % Lung erosion
                bt3 = 9; % End dilation
                bt4 = 3; % End erosion
                slc = round(2*self.dims(3)/3);
            case 3 % Clin. Exp.
                sm = 3;
                b = self.scaleB(vec);
                Ti = 607 + b;
                Tt = 125 + b;
                bt1 = 8; % Trachea dilation
                bt2 = 8; % Lung erosion
                bt3 = 9; % End dilation
                bt4 = 3; % End erosion
                slc = round(2*self.dims(3)/3);
            case 4 % Mouse Insp.
                % Low-res
%                 sm = 0;
%                 b = self.scaleB(vec);
%                 Ti = 500+b;
%                 Tt = 200+b;
%                 bt1 = 1; % Trachea dilation
%                 bt2 = 4; % Lung erosion
%                 bt3 = 5; % End dilation
%                 bt4 = 1; % End erosion
%                 slc = round(2*self.dims(3)/3);
                % High-res
                sm = 4;
                b = self.scaleB(vec);
                Ti = 500+b;
                Tt = 150+b;
                bt1 = 3; % Trachea dilation
                bt2 = 4; % Lung erosion
                bt3 = 5; % End dilation
                bt4 = 2; % End erosion
                slc = round(2*self.dims(3)/3);
            case 5 % Mouse Exp.
                sm = 4;
                b = self.scaleB(vec);
                Ti = 700;
                Tt = 300;
                bt1 = 8; % Trachea dilation
                bt2 = 6; % Lung erosion
                bt3 = 5; % End dilation
                bt4 = 2; % End erosion
                slc = round(0.7*self.dims(3));
        end
        dopt = ~isempty(slc);

        disp('.')
        disp('.')
        disp('Starting Lung Segmentation ...')
        tic

        if sm
            tmat = filtGaussSep(self.mat(:,:,:,vec),sm);
        else
            tmat = self.mat;
        end
        % Make voxels outside FOV have air HU value
        tmat(tmat<b) = b;

        % First, find optimal lung threshold
        if isempty(Ti)
            To = 0; ct = 0; maxct = 100;
            disp('  Finding threshold ...')
            while (abs(Ti-To)>0.5) && (ct<=maxct)
                To = Ti;
                ub = mean(tmat(tmat(:)>=To));
                un = mean(tmat(tmat(:)<To));
                Ti = (ub + un)/2;
                ct = ct+1;
                disp(['    Count: ' num2str(ct) ' ; Ti: ' num2str(round(Ti))])
            end
        end
        if isempty(Tt)
            % then find auto air threshold
            To = 0; ct = 0;
            tmat2 = tmat(tmat<Ti);
            disp('  Finding trach. threshold ...')
            while (abs(Tt-To)>0.5) && (ct<=maxct)
                To = Tt;
                ub = mean(tmat2(tmat2(:)>=To));
                un = mean(tmat2(tmat2(:)<To));
                Tt = (0.7*ub + un)/2;
                ct = ct+1;
                disp(['    Count: ' num2str(ct) ' ; Tt: ' num2str(Tt)])
            end
        end

        % Find image ends
        tmask = false(self.dims(1:3));
        tmask([1,end],:,:) = true;
        tmask(:,[1,end],:) = true;
        cind = find(tmask);
        tmask = self.mask.mat;
        vind = find(self.mask.mat);

        % Create initial filtered masks for lung and trachea
        if dopt
            th = tic;
        end
        % Trachea mask
        if pcheck,self.mask.merge('Replace',tmat<Tt);disp('Line 190'),pause,end
        lmask = filtGaussSep(double(tmat>=Tt),2) <= 0.1; % erode to get rid of small regions
        if pcheck,self.mask.merge('Replace',lmask);disp('Line 192'),pause,end
        lmask = filtGaussSep(lmask,bt1) > 0.1; % dilate to fill trachea
        
        % Current VOI marks bottom of bronchi to cut
        %ind = find(any(any(tmask,1),2),1);
        %lmask(:,:,1:ind-1) = false;
        
        cc = bwconncomp(lmask);
        lmask = false(size(lmask));
        for i = 1:cc.NumObjects
            if any(ismember(vind,cc.PixelIdxList{i}))
                lmask(cc.PixelIdxList{i}) = true;
            end
        end
        
        lmask = filtGaussSep(double(lmask),bt1) > 0.01; % dilate
        if pcheck,self.mask.merge('Replace',lmask);pause,end
        if dopt
            hf=figure('Position',[40 40 600 800]);
            subplot(3,2,1),
            imshow(lmask(:,:,slc));
            title('Trachea')
            disp(['  Filterd trachea mask (' num2str(toc(th)) ' seconds)...'])
            pause(0.1)
            th = tic;
        end
        % Lung mask - w/o trachea
        lmask = filtGaussSep(1-double((tmat<Ti) & ~lmask),bt2) < 0.1; % erode
        if pcheck,self.mask.merge('Replace',lmask);pause,end
        if dopt
            figure(hf),subplot(3,2,2),
            imshow(lmask(:,:,slc));
            title('Lungs')
            disp(['  Filtered lung mask (' num2str(toc(th)) ' seconds)...'])
            pause(0.1)
            th = tic;
        end
        %clear tmat

        % Find connected regions
        cc = bwconncomp(lmask);
        if dopt
            disp(['  Found ' num2str(cc.NumObjects) ' regions (' num2str(toc(th)) ' seconds)...'])
        end

        % Discard regions that are too small and largest (outside of body)
        n = numel(lmask);
        nvox = cellfun(@numel,cc.PixelIdxList);
        ind = cellfun(@(x)any(ismember(x,cind)),cc.PixelIdxList);
        idx = find(~ind & ((nvox/n)>0.005));
        disp(['    Using ' num2str(length(idx)) ' regions ...'])
        lmask2 = false(size(lmask)); lmask = lmask2;
        for i = 1:length(idx)
            disp(['    Region ' num2str(i) ': ' num2str(nvox(idx(i))) ' voxels (',...
                num2str(nvox(idx(i))*prod(self.voxsz)/10^6) ' L)'])
            if dopt
                th2 = tic;
            end
            lmask2 = false(size(lmask));
            lmask2(cc.PixelIdxList{idx(i)}) = true;
            lmask2 = filtGaussSep(double(lmask2),bt3) > 0.1; % dilate
            if dopt
                figure(hf),subplot(3,2,2*i+1),
                imshow(lmask2(:,:,slc));
                title('Dilated')
                disp(['      Dilated (' num2str(toc(th2)) ' seconds)...'])
                pause(0.1)
                th2 = tic;
            end
            cc2 = bwconncomp(~lmask2);
            nvox2 = cellfun(@numel,cc2.PixelIdxList);
            lmask2 = true(size(lmask));
            lmask2(cc2.PixelIdxList{nvox2==max(nvox2)}) = false;
            if dopt
                disp(['      Found outside (' num2str(toc(th2)) ' seconds)...'])
                th2 = tic;
            end
            lmask2 = filtGaussSep(1-double(lmask2),bt4) < 0.1; % erode
            if dopt
                figure(hf),subplot(3,2,2*i+2),
                imshow(lmask2(:,:,slc));
                title('Eroded')
                disp(['      Eroded (' num2str(toc(th2)) ' seconds)...'])
                pause(0.1)
            end
            lmask = (lmask | lmask2);
        end
        disp(['Ended after ' num2str(toc) ' seconds'])

        self.mask.merge('Replace',lmask);
    end
else
    errordlg('Need VOI to select trachea / major bronchi!')
end