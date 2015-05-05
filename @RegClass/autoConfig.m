% RegClass function
function autoConfig(self, ~, ~)
% Function for automatically determining parameters for coregistration
% based on image dimensions

if self.cmiObj(1).img.check && self.cmiObj(2).img.check
    
    d1 = self.cmiObj(1).img.dims(1:3);
    d2 = self.cmiObj(2).img.dims(1:3);
    fov1 = d1.*self.cmiObj(1).img.voxsz;
    fov2 = d2.*self.cmiObj(2).img.voxsz;
    
    % Fixed Image Pyramid
    f = round(d1/64);
    
    % Moving Image Pyramid
    m = round(d2/64);
    
    C = {'NumberOfResolutions',nres,...
        'FixedImagePyramidSchedule',f,...
        'MovingImagePyramidSchedule',m};
    
    % SP_a
    SPa = 2*sqrt(sum(d2.^2))^(1/3);
    
    % Grid
    if strcmp(self.elxObj.getPar(self.ind,'Transform'),'BSplineTransform')
        SPa = 
        gsp
        if max([fov1,fov2])>100
            be = 50; % Human data
        else
            be = 0.1; % Preclinical
        end
        C = [C,{'GridSpacingSchedule',gsp,'TransformBendingEnergy',be}];
    end
    
    % Image Sampler
    if self.cmiObj(1).img.mask.check
        smpler = 'RandomSparseMask';
    else
        smpler = 'RandomCoordinate';
    end
    C = [C,{'ImageSampler',smpler}];
    
    % Set the parameters
    self.setElxPar(self.ind,C{:});
end