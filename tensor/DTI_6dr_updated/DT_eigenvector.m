function EVector=DT_eigenvector(EV,D)
% [EV,IVmask]=DT_eigenvalue_B(D); % Calculate eigenvalues
% % Calculate eigenvalues, eigenvectors, and other terms
% Ref: Hasan et al., JMR (2001) 152:41-47.

warning off MATLAB:divideByZero

% EV=DT_eigenvalue_B(D);

A=zeros(size(D,1),size(D,2),size(D,3));
B=zeros(size(D,1),size(D,2),size(D,3));
C=zeros(size(D,1),size(D,2),size(D,3));

EVector1=zeros(size(D,1),size(D,2),size(D,3),1,3);
EVector=zeros(size(D,1),size(D,2),size(D,3),1,3);
EVectn=zeros(size(D,1),size(D,2),size(D,3),1);

'Calculating Eigenvectors: DT_eigenvector'

for i=1:3
    A(:,:,:)=D(:,:,:,1)-EV(:,:,:,i);
    B(:,:,:)=D(:,:,:,4)-EV(:,:,:,i);
    C(:,:,:)=D(:,:,:,6)-EV(:,:,:,i);
    
    EVector1(:,:,:,1,1)=(D(:,:,:,2).*D(:,:,:,5)-B(:,:,:).*D(:,:,:,3)).*(D(:,:,:,3).*D(:,:,:,5)-C(:,:,:).*D(:,:,:,2)); % eix matlab=columns [0 1 0]
    EVector1(:,:,:,1,2)=(D(:,:,:,3).*D(:,:,:,5)-C(:,:,:).*D(:,:,:,2)).*(D(:,:,:,3).*D(:,:,:,2)-A(:,:,:).*D(:,:,:,5)); % eiy matlab=rows [1 0 0]
    EVector1(:,:,:,1,3)=(D(:,:,:,2).*D(:,:,:,5)-B(:,:,:).*D(:,:,:,3)).*(D(:,:,:,3).*D(:,:,:,2)-A(:,:,:).*D(:,:,:,5)); % eiz matlab=slices [0 0 1]
    
    % When using the masks, the artifacts are seen in the quiver plot.
    EVectn(:,:,:,1)=(sqrt(EVector1(:,:,:,1,1).^2+EVector1(:,:,:,1,2).^2+EVector1(:,:,:,1,3).^2));
    
    EVector(:,:,:,1,1)=(EVector1(:,:,:,1,1)./EVectn(:,:,:,1)); % X component
    EVector(:,:,:,1,2)=(EVector1(:,:,:,1,2)./EVectn(:,:,:,1)); % Y component
    EVector(:,:,:,1,3)=(EVector1(:,:,:,1,3)./EVectn(:,:,:,1)); % Z component
    
    EVector(:,:,:,:,1)=-EVector(:,:,:,:,1);  % Note:  Negative sign is so calculation of EV conforms with MATLAB
    EVector(:,:,:,:,2)=EVector(:,:,:,:,2);
    EVector(:,:,:,:,3)=EVector(:,:,:,:,3);
    
end
