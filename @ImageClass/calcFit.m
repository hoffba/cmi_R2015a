% ImageClass function
% Perform exponential/perfusion/diffusion analysis on image data
function calcFit(self)
if self.check && self.mask.check
    stat = true;
    nx = self.dims(4); % number of data points on curve
    % Check to see if model XData needs to be updated
    if length(self.model.x)~=self.dims(4)
        x = cellfun(@str2double,self.labels);
        if any(isnan(x)) % ask for input
            stat = false;
            while ~stat % loop until either cancelled or valid input
                stat = self.setFitX;
            end
        else % use labels to set XData
            stat = self.model.setXData(x);
        end
    end
    if stat
        if self.voxfit % Perform a voxel-wise fit
        % Extract VOI voxel data
            inds = find(self.mask.mat(:));
            npts = length(inds);
            nvol = prod(self.dims(1:3));
            ntot = prod(self.dims);
        % Perform first voxel fit
            npar = length(self.model.getPar0);
            nvout = npar+2;
            pdata = zeros(npts,npar);
            gof = zeros(npts,1);
        % Loop over the voxels
            hw = waitbar(0);
            wstrs = {'Voxel # ',[' / ',num2str(npts)]};
            for i = 1:npts
                waitbar(i/npts,hw,[wstrs{1},num2str(i),wstrs{2}])
                self.model.setYData(self.mat(inds(i):nvol:ntot));
                try
                    pdata(i,:) = self.model.fitData;
                catch err
                    rethrow(err);
                end
                gof(i) = self.model.calcGoF(pdata(i,:));
            end
            delete(hw);
        % Replace current image with fit maps
            bgLabel = {};
            plabels = self.model.getDefs.labels;
            switch self.model.getModType
                case 'diff'
                    % Keep lowest b-value image for anatomical ref
                    bgLabel = self.labels(1);
                    bgvec = 1;
                    nvf = 1;
                case 'perf' % DCE Model
                    % Keep a T1 post-Gd image for anatomical ref
                    bgLabel = {'T1post'};
                    bgvec = nx;
                    plabels = [{'AUC'},plabels];
                    nvout = nvout + 1; % for AUC
                    nvf = 2;
                    bpts = self.model.dceOpts.nbase;
                    auc = sum(self.mat,4) ./ ...
                        (bpts*mean(self.mat(:,:,:,1:bpts),4));
                    auc(isinf(auc)) = 0;
                    self.mat(:,:,:,2) = auc;
            end
            % First image (for background / anatomy)
            self.mat(:,:,:,1) = self.mat(:,:,:,bgvec);
            % Delete unnecessary 4D images
            self.imgDelete((nvout+1):nx);
            % Next, set parameter maps
            for i = 1:npar
                tmap = 0*self.mask.mat;
                tmap(self.mask.mat) = pdata(:,i);
                self.mat(:,:,:,i+nvf) = tmap;
            end
            % Last, set GoF map
            tmap = 0*self.mask.mat;
            tmap(self.mask.mat) = gof;
            self.mat(:,:,:,nvout) = tmap;
            % Update image-related parameters
            self.labels = [bgLabel,plabels,{'GoF'}];
            self.thresh = 100000*([-1;1]*ones(1,nvout))';
            self.scaleM = ones(1,nvout);
            self.scaleB = zeros(1,nvout);
        else % Perform a volume-mean fit
            ydata = zeros(1,nx);
            for i = 1:nx
                td = self.mat(:,:,:,i);
                ydata(i) = mean(td(self.mask.mat));
            end
            self.model.setYData(ydata);
            pout = self.model.fitData;
            pLabels = self.model.getDefs.labels;
            np = length(pout);
            txt = cell(np,1);
            for i = 1:np
                txt{i} = [pLabels{i} '  =  ' num2str(pout(i))];
            end
            helpdlg(txt);
        end
    end
end