% ImageClass function
% Calculate PRM
function [labels,vals] = calcPRM(self,vec)
labels = {}; vals = [];
% Make sure image and mask are available
if self.check && self.mask.check
    if (nargin == 2) && isnumeric(vec)
        
        % extract mask coordinates
        mskind = find(self.mask.mat);
        np = numel(mskind);
        
        % determine what images need to be passed in
        tvec = self.prm.dvec;
        tvec(tvec==0) = vec;
        tvec = unique(tvec);
        nv = length(tvec);
        
%         if nv>1
            % extract image values
            vals = zeros(np,nv);
            for i = 1:nv
                vals(:,i) = self.mat(mskind + (tvec(i)-1)*(prod(self.dims(1:3))));
            end

            % Finally, calculate the PRM
            [labels,vals] = self.prm.calcPRM(vals,tvec,vec,self.labels(tvec),self.mask.mat);
%         end
    end
end