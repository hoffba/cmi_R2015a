% ImageClass function
% Calculate PRM
function [labels,vals] = calcPRM(self,vec)
labels = {}; vals = [];
% Make sure image and mask are available
if self.check && self.mask.check
    if (nargin == 2) && isnumeric(vec)
        
        % determine what images need to be passed in
        tvec = self.prm.dvec;
        tvec(tvec==0) = vec;
        tvec = unique(tvec);
        
        if self.prm.normchk
            self.imgFilt(tvec,'median',[3,3]);
        end
        
        % extract image values
        vals = self.getMaskVals(tvec);

        % Finally, calculate the PRM
        [labels,vals] = self.prm.calcPRM(vals,tvec,vec,self.labels(tvec),self.mask.mat);
    end
end