function [D ADC Trace EValue EVector FA]=DTI_analysis(data);
%
%
% to get DTI results:
% 1. load fdf DTI data using CMI program.
% 2. type data=cmiObj0.img.mat at prompt
% 3. run this program at prompt ">>[D ADC Trace EValue EVector
% FA]=DTI_analysis(data);"


%%
param{1}=1000; % b values
param{2}=6; % diffusion scheme hardcoded for Agilent scheme


%% Diffusion Scheme Agilent scheme
% param{3}=[[0 0 0];[1 1 0];[1 0 1];[0 1 1];[-1 1 0];[-1 0 1];[0 -1 1]];


%%
bV=param{1};


%% Calculate Normalized Signal Intesity
IM=zeros(size(data(:,:,:,1:6)));
for i=1:size(data,4)-1
    IM(:,:,:,i)=log(data(:,:,:,i+1)./data(:,:,:,1));
end

IM=IM(:,:,:,[2,5,3,6,1,4]);
%% Calculate Algebraic solution of Diffusion Tensor
% Basser & Pierpaoli, MRM, (1996) 39:928-934
D=zeros(size(IM));
D(:,:,:,1)=((-IM(:,:,:,1)-IM(:,:,:,2)+IM(:,:,:,3)+IM(:,:,:,4)-IM(:,:,:,5)-IM(:,:,:,6))./(2*bV)); % Dxx
D(:,:,:,2)=((-IM(:,:,:,5)+IM(:,:,:,6))./(2*bV)); % Dxy
D(:,:,:,3)=((-IM(:,:,:,1)+IM(:,:,:,2))./(2*bV)); % Dxz
D(:,:,:,4)=((IM(:,:,:,1)+IM(:,:,:,2)-IM(:,:,:,3)-IM(:,:,:,4)-IM(:,:,:,5)-IM(:,:,:,6))./(2*bV)); % Dyy
D(:,:,:,5)=((-IM(:,:,:,3)+IM(:,:,:,4))./(2*bV)); % Dyz
D(:,:,:,6)=((-IM(:,:,:,1)-IM(:,:,:,2)-IM(:,:,:,3)-IM(:,:,:,4)+IM(:,:,:,5)+IM(:,:,:,6))./(2*bV)); % Dzz

%% Eigenvalues and Eigenvectors
[EValue,ADC,Trace,FA]=DT_eigenvalue(D);
EVector=DT_eigenvector(EValue,D); % not sure about orientation of eigenvectors
orien1=repmat(FA,[1,1,1,3]).*abs(squeeze(EVector(:,:,:,1,[2 1 3]))); % FA*EVector=orientation map resorted
%% Temp figures
% figure(99);imshow(squeeze(orien1(:,:,8,:));
% figure(100);imagesc(FA(:,:,17));
[d1 d2 d3]=size(ADC);
if d1==128&&d2==128
    ADC=imresize(ADC,2);FA=imresize(FA,2);orien1=imresize(orien1,2);
end

saveFLD('DTI_proc',cat(4,ADC*10^6,FA*10^3),[{'ADCe6'},{'FAe3'}],size(ADC)); % save data as fld
saveFLD('DT_orien',orien1*10^3,{'Orien1e3'},size(orien1));
assignin('base','Evector',EVector);
assignin('base','FA',FA);
assignin('base','orien1',orien1);

    
