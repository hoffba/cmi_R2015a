% ImageClass function
% perform perfusion analysis on MRI data
function analyzePerf(self)
if self.check && self.mask.check
    nx = self.dims(4); % number of data points on curve
    if self.voxfit % Perform a voxel-wise fit
    % Extract VOI voxel data
        inds = find(self.mask.mat(:));
        npts = length(inds);
        nvol = prod(self.dims(1:3));
        ntot = prod(self.dims);
    % Perform first voxel fit
        self.model.setYData(self.mat(inds(1):nvol:ntot));
        [pout,pLabels] = self.model.fitData;
        npar = length(pout);
        pdata = zeros(npts,npar);
        pdata(1,:) = pout;
    % Loop over the rest of the voxels
        for i = 2:npts
            self.model.setYData(self.mat(inds(i):nvol:ntot));
            [pdata(i,:),~] = self.model.fitData;
        end
    % Replace current image with fit maps
        mt = self.model.getModelType;
        bgLabel = {};
        switch mt
            case mt==1
                % Keep lowest b-value image for anatomical ref
                bgLabel = self.labels(1);
                bgvec = 1;
                npar = npar+1;
            case mt==2 % DCE Model
                % Keep a T1 post-Gd image for anatomical ref
                bgLabel = {'T1post'};
                bgvec = nx;
                npar = npar+1;
        end
        self.mat(:,:,:,1) = self.mat(:,:,:,bgvec);
        self.imgDelete((npar+1):nx);
        self.labels = [bgLabel pLabels];
        for i = 2:npar
            tmap = 0*self.mask.mat;
            tmap(self.mask.mat) = pdata(:,i-1);
            self.mat(:,:,:,i) = tmap;
        end
        self.thresh = 100000*([-1;1]*ones(1,npar))';
        self.scaleM = ones(1,npar);
        self.scaleB = zeros(1,npar);
    else % Perform a volume-mean fit
        nt = self.dims(4);
        ydata = zeros(1,nt);
        for i = 1:nt
            td = self.mat(:,:,:,i);
            ydata(i) = mean(td(self.mask.mat));
        end
        self.model.setYData(ydata);
        [pout,pLabels] = self.model.fitData;
        np = length(pout);
        txt = cell(np,1);
        for i = 1:np
            txt{i} = [pLabels{i} '  =  ' num2str(pout(i))];
        end
        helpdlg(txt);
    end
end