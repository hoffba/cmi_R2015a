function I = tensorShape(E,str)

if (nargin==2) && (size(E,4)==3)
    if ischar(str)
        str = {str};
    end
    if iscellstr(str)
        d = size(E);
        npar = length(str);
        
        if any(ismember({'Det','VR','VF'},str))
            D = prod(E,4);
        end
        if any(ismember({'CL','CP','CS','CA','RA','FA','VR','VF'},str))
            M = mean(E,4);
        end

        I = nan([d(1:3),npar]);
        for i = 1:npar
            switch str{i}
                case 'Det'
                    tI = D;
                case 'ADI' % Anisotropic deformation index
                    tI = sqrt(((E(:,:,:,1)-E(:,:,:,2))./E(:,:,:,2)).^2 ...
                            + ((E(:,:,:,2)-E(:,:,:,3))./E(:,:,:,3)).^2);
                    tI(tI>100) = 100;
                case 'SRI' % Slab-rod index
                    tI = atan( E(:,:,:,3) .* (E(:,:,:,1)-E(:,:,:,2)) ...
                            ./ E(:,:,:,2) ./ (E(:,:,:,2)-E(:,:,:,3)) )/pi*2;
                case 'CL' % Linear coefficient
                    tI = (E(:,:,:,1) - E(:,:,:,2)) ./ (3*M);
                case 'CP' % Planar coefficient
                    tI = 2*(E(:,:,:,2)-E(:,:,:,3)) ./ (3*M);
                case 'CS' % Spherical coefficient
                    tI = E(:,:,:,3)./M;
                case 'CA' % Anisotropy coefficient
                    tI = 1 - E(:,:,:,3)./M;
                case 'RA' % Relative anisotropy
                    tI = sqrt(sum((E-repmat(M,1,1,1,3)).^2,4))./(sqrt(6)*M);
                case 'FA' % Fractional anisotropy
                    tI = sqrt(3*sum((E-repmat(M,1,1,1,3)).^2,4)./(2*sum(E.^2,4)));
                case 'VR' % Volume ratio
                    tI = D./(M.^3);
                case 'VF' % Volume fraction
                    tI = ((M.^3) - D)./(M.^3);
            end
            I(:,:,:,i) = tI;
        end
    end
end