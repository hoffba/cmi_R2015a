% ImageClass function
% Re-orient image based on Eigenvector analysis of current VOI (for CT bone)
function eigOrient(self)
    dorder = [3 2 1];
    %% VOI -> indices
        [X,Y,Z] = self.getImageCoords;
        pcain = [X(self.mask.mat),Y(self.mask.mat),Z(self.mask.mat)];
    %% PCA analysis
        p = size(pcain,2);
        [score,sigma,coeff] = svd(pcain,0);
        sigma = diag(sigma);
        score = bsxfun(@times,score,sigma');
        [~,maxind] = max(abs(coeff),[],1); % Enforce a sign convention, largest component - positive
        d = size(coeff,2);
        colsign = sign(coeff(maxind + (0:p:(d-1)*p)));
        score = bsxfun(@times,score,colsign);
    %% Map to new space
        T = self.calcTxF(pcain,score(:,dorder));
        self.affTxF(T,'linear');
        self.mask.affTxF(T,'linear');
%end