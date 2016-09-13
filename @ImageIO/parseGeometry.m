% ImageIO method
function geo = parseGeometry(self,info,itype)
% Input raw image info from loader functions
% Returns structure containing CMI-required information

geo = self.definfo;
switch itype
    case 1 % Analyze
    case 2 % MetaIO
        
    case 3 % AVS
    case 4 % FDF
    case 5 % VFF
    case 6 % DICOM
    case 7 % NifTi
    case 8 % Bruker
    case 9 % MRSolutions
    case 10 % Mask
    case 11 % TIFF
    case 12 % JPEG
    case 13 % Matlab
    case 14 % FID
end