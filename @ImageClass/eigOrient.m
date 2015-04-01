% ImageClass function
% Re-orient image based on Eigenvector analysis of current VOI (for CT bone)
function eigOrient(self)
    dorder = [3 2 1];
    %% VOI -> indices
        [y,x,z] = ind2sub(self.dims(1:3),find(self.mask.mat));
        pcain = [x,y,z];
    %% PCA analysis
        pcamean = mean(pcain,1);
        p = size(pcain,2);
        score = bsxfun(@minus,pcain,pcamean); % Center the image coordinates on the origin
        [score,sigma,coeff] = svd(score,0);
        sigma = diag(sigma);
        score = bsxfun(@times,score,sigma');
        [~,maxind] = max(abs(coeff),[],1); % Enforce a sign convention, largest component - positive
        d = size(coeff,2);
        colsign = sign(coeff(maxind + (0:p:(d-1)*p)));
        score = bsxfun(@times,score,colsign);
        %score = bsxfun(@plus,score(:,dorder),pcamean(dorder));
        score = bsxfun(@plus,score(:,dorder),round(self.dims(1:3)/2));
    %% Map to new space
        T = self.calcTxF(pcain,score);
        self.affTxF(T,'linear');
        self.mask.affTxF(T,'linear');
%end