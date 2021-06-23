
function img = grid2img(MF,ind,mask,nd,stat)
% Interpolate grid data to original dimensions
% Inputs:
%   MF = [nth,nmf,np] vector of values
%   ind = [npx1] vector matrix indices
%   mask = binary mask 3D matrix

img = [];
if nargin>2
    
    [nth,nmf,np] = size(MF);
    
    % Validate inputs:
    if np~=length(ind)
        error('Values and indices must be same size.');
    end
    if ~all((ind>0)&(ind<=numel(mask)))
        error('Invalid indices.');
    end
    if nargin<4
        nd = 3;
    end
    gchk = true;
    ssz = get(0,'ScreenSize');
    if all(ssz==1)
        gchk = false; % This indicates batch mode
    end
    
    if (nargin<5) || (stat==0)
    % Start batch:
    
        img = batch(@grid2img,1,{MF,ind,mask,nd,1});
        
    else
    % Process data:
    
        d = size(mask);
        img = zeros([d,nmf,nth]);
        switch nd % Either do 2D or 3D interpolation
            case 2 % mainly for slice-gapped data
                disp('Calculating meshgrid ...')
                [Y,X,Z] = ind2sub(d,ind);
                % Loop only over slices with VOI:
                zi = intersect( find(any(any(mask,1),2)) , Z );
                nz = length(zi);
                [Xq,Yq] = meshgrid(1:d(2),1:d(1));
                if gchk
                    hw = waitbar(0,'Interpolating slice ...');
                end
                for ith = 1:nth
                    for imf = 1:nmf
                        timg = zeros(d);
                        for i = 1:nz
                            if gchk
                                waitbar(i/nz,hw,sprintf('Interpolating slice %u',zi(i)));
                            else
                                fprintf('Intepolating: th%u, mf%u, slc%u\n',ith,imf,i);
                            end
                            
                            ii = (Z==zi(i));
                            if nnz(ii)>2
                                F = scatteredInterpolant(X(ii),Y(ii),squeeze(MF(ith,imf,ii)),'linear','none');
                                timg(:,:,zi(i)) = F(Xq,Yq);
                            end
                        end
                        timg(~mask | isnan(timg)) = 0;
                        img(:,:,:,imf,ith) = timg;
                    end
                end
                if gchk, delete(hw); end
            case 3
                % Generate interpolation model:
                disp('Calculating meshgrid ...')
                [Xq,Yq,Zq] = meshgrid(1:d(2),1:d(1),1:d(3));
                
                [Y,X,Z] = ind2sub(d,ind);
                for ith = 1:nth
                    for imf = 1:nmf
                        fprintf('Intepolating: th%u, mf%u\n',ith,imf);
                        F = scatteredInterpolant(X,Y,Z,squeeze(MF(ith,imf,:)),'linear','none');
                        % Calculate full grid
                        timg = F(Xq,Yq,Zq);
                        timg(~mask | isnan(timg)) = 0;
                        img(:,:,:,imf,ith) = timg;
                    end
                end
            otherwise
                error('Invalid number of dimensions. Must be 2 or 3');
        end
    end
end





