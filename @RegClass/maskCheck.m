function [ check, mask ] = maskCheck(self)
%Checks if mask is loaded on either of the two images and returns a true
%value if such a mask is loaded and the value of the mask itself.
%check if mask is loaded or not
if ~self.cmiObj(2).img.mask.check && ~self.cmiObj(1).img.mask.check
    yn = questdlg('Mask is not loaded on either object, would you like to continue anyways?', '404 mask not found', 'Yes', 'No');
elseif self.cmiObj(2).img.mask.check && ~self.cmiObj(1).img.mask.check
    mask = self.cmiObj(2).img.mask.mat;
    yn = 'Yes';
elseif self.cmiObj(1).img.mask.check && ~self.cmiObj(2).img.mask.check
    mask = self.cmiObj(1).img.mask.mat;
    yn = 'Yes';
elseif self.cmiObj(1).img.mask.check && self.cmiObj(2).img.mask.check
    mskChoice = questdlg('Mask is loaded on both objects, which one would you like to use?', 'Choose a mask', 'Reference', 'Homologous');
    switch(mskChoice)
        case 'Reference'
            mask = self.cmiObj(1).img.mask.mat;
            yn = 'Yes';
        case 'Homologous'
            mask = self.cmiObj(2).img.mask.mat;
            yn = 'Yes';
    end
end
check = strcmp(yn, 'Yes');

end

