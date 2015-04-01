function stat = saveAffineTxformPars(pars,info,fname)
% Save Elastix transform parameter file
% - Assumes 3D image transform
% Inputs:
%   pars = vector of elastix transform parameters
%       Length(pars) determines transform type:
%       Translation (3) : [tx,ty,tz]
%       Euler       (6) : [rx,ry,rz,tx,ty,tz]   * r is the angle (radians)
%       Similarity  (7) : [qx,qy,qz,tx,ty,tz,s] * q is the rotational versor
%       Affine      (12): [a11,a12,a13,a21,a22,a23,a31,a32,a33,tx,ty,tz]
%   info = struct: {dims , voxsz} (of reference geometry)
%   fname = full file path of desired text file

stat = false;
if (nargin>1)
    % Determine transformation type:
    % Assume 3D image
    np = length(pars);
    switch np
        case 3 % Translation
            Ttype = 'TranslationTransform';
        case 6 % Euler
            Ttype = 'EulerTransform';
        case 7 % Similarity
            Ttype = 'SimilarityTransform';
        case 12 % Affine
            Ttype = 'AffineTransform';
        otherwise
            error(['Number of parameters (',num2str(np),...
                    ') does not match any types.']);
    end
    % Write Transformation file:
    if (nargin<3)
        uiputfile();
    end
    fid = fopen(fname,'w');
    if fid==-1
        error('Could not open file to write.');
    else
        fprintf(fid,'(Transform "%s")\n',Ttype);
        fprintf(fid,'(NumberOfParameters %u)\n',length(pars));
        fprintf(fid,'(TransformParameters%s)\n',sprintf(' %.6f',pars));
        fprintf(fid,'(InitialTransformParametersFileName "NoInitialTransform")\n');
        fprintf(fid,'(HowToCombineTransforms "Compose")\n');
        
        fprintf(fid,'\n// Image specific\n');
        fprintf(fid,'(FixedImageDimension 3)\n');
        fprintf(fid,'(MovingImageDimension 3)\n');
        fprintf(fid,'(FixedInternalImagePixelType "float")\n');
        fprintf(fid,'(MovingInternalImagePixelType "float")\n');
        fprintf(fid,'(Size%s)\n',sprintf(' %u',info.dims));
        fprintf(fid,'(Index%s)\n',sprintf(' %u',zeros(1,3)));
        fprintf(fid,'(Spacing%s)\n',sprintf(' %.10f',info.voxsz));
        fprintf(fid,'(Origin%s)\n',sprintf(' %.10f',-(info.dims-1).*info.voxsz/2));
        fprintf(fid,'(Direction%s)\n',sprintf(' %u',[1,0,0,0,1,0,0,0,1]));
        fprintf(fid,'(UseDirectionCosines "%s")\n','false');
        
        if ismember(np,[6,7,12])
            fprintf(fid,'\n// %s specific\n',Ttype);
            fprintf(fid,'(CenterOfRotationPoint%s)\n',...
                    sprintf(' %u',zeros(1,3)));
        end
        if np==6
            fprintf(fid,'(ComputeZYX "false")\n');
        end
        
        fprintf(fid,'\n// ResampleInterpolator specific\n');
        fprintf(fid,'(ResampleInterpolator "%s")\n','FinalBSplineInterpolator');
        fprintf(fid,'(FinalBSplineInterpolationOrder %u)\n',3);
        
        fprintf(fid,'\n// Resampler specific\n');
        fprintf(fid,'(Resampler "DefaultResampler")\n');
        fprintf(fid,'(DefaultPixelValue %.6f)\n',0);
        fprintf(fid,'(ResultImageFormat "mhd")\n');
        fprintf(fid,'(ResultImagePixelType "float")\n');
        fprintf(fid,'(CompressResultImage "false")\n');
        stat = fclose(fid);
        stat = ~stat;
    end
end