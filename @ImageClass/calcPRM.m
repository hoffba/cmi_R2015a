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
        nvec = length(tvec);
        
        % Optional pre-filter of image, then grab masked values:
        if self.prm.filtchk
            timg = ImageClass;
            timg.setMat(self.mat(:,:,:,tvec),...
                cellfun(@num2str,num2cell(tvec),'UniformOutput',false),...
                self.voxsz.*self.dims(1:3),self.orient,self.name); %CJG added 20210812
            timg.mask.merge('replace',self.mask.mat);
            timg.imgFilt(1:nvec,self.prm.filttype,{self.prm.filtstr});
            vals = timg.getMaskVals;
            delete(timg);
        else
            vals = self.getMaskVals(tvec);
        end

        % Finally, calculate the PRM
        [labels,vals] = self.prm.calcPRM(vals,tvec,vec,self.labels(tvec),self.mask.mat);
    end
end