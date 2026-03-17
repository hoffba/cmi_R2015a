function V = vesselSeg_subj(ct, lobes)

% Scale lung image to between 0 and 1
    I = (single(ct) + 1000)/1000;
    I(~lobes | I>1) = 1;
    I(I<0 | isnan(I)) = 0;
    
    %MTHT3D Parameters
    no = 12;
    beta = 70; 
    alpha = 0.5; 
    c = 15; %parameters for the Vesselness
    scale = [0.5,1,1.5,3:2:13];

    V = vessel_MTHT3D(I,lobes,scale,no,beta,alpha,c);

end

function V = vessel_MTHT3D(I,lobes,s,no,beta,alpha,c)
%%
%	INPUT:
%       I     - 3D input image (image normalized to between 0,1)
%		s     - Scale
% 		beta  - Parameter for vesselness
%       alpha - Parameter for vesselness
%		c     - Parameter for vesselness
%
%   OUTPUT:
%       V - Vesselness Enhanced Image
%
%   AUTHOR:
%       Sundaresh Ram
%       Benjamin Hoff
% 		Galban Lab
% 		Dept. Of Radiology & Dept. of BME
% 		University of Michigan
% 		e-mail : sundarer@umich.edu
%
%   VERSION:
%       0.1 - 02/21/2020 First implementation
%       1.0 - 03/12/2026 Adjusted for implementation in the MiTAP pipeline
%
    d = size(I,1:3);
    idx = find(lobes);
    np = length(idx);
    Vmax = single(0);

    % Point distribution on the sphere of unit radius
    [~,orients,~] = SurfacesSpiralPoints3D(no);
    for i = 1:length(s) % Loop over scales

        fprintf(['Scale: ' num2str(s(i)),'\n']);
        T = single(0);
        for j = 1:size(orients,1) % Loop over orientations
            se = Line3D(s(i),orients(j,:),1,1,1); %3D line structuring element
            im_temp = imopen(I,se);

            % Tensor
            nnT = reshape(kron(orients(j,:),orients(j,:)')',1,[]);
            T = T + (I(idx)-im_temp(idx)) .* nnT;
        end

        % Eigen matrix
        T = permute(reshape(T,[np,3,3]),[3,2,1]);
        T = sort(eig3(T),1); % [3 x np] : [L1;L2;L3] Eigenvalues

        % Vesselness
        Vmax = max(Vesselness(beta,alpha,c,T),Vmax);
    end

    V = zeros(d);
    V(idx) = Vmax;
end

function  V = Vesselness(beta,alpha,c,L)
    Ralpha = abs(L(2,:))./abs(L(3,:));              % 1 -> plate;   0 -> line
    Rbeta = abs(L(1,:))./sqrt(abs(L(2,:).*L(3,:))); % 1 -> blob;    0 -> line
    Ralpha(isnan(Ralpha)) = 0;
    Rbeta(isnan(Rbeta)) = 0;
    S = sqrt(L(1,:).^2 + L(2,:).^2 + L(3,:).^2);     
    V =  ( 1 - exp(-(Ralpha.^2)/(2*alpha^2)) ) .* exp(-(Rbeta.^2)/(2*beta^2)) .*  ( 1 - exp(-(S.^2)/(2*c^2)) );
    V(any(L<=0,1)) = 0;     
    V = Normalize(V);
end  