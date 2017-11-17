function [EValue,ADC,Trace,FA]=DT_eigenvalue(D)
'DT_eigenvalue'

I=DT_invariant(D); % Calculate invariants
% Calculate eigenvalues, eigenvectors, and other terms

% 'Calculate eigenvalues'
v1=(I(:,:,:,1)/3).^2-I(:,:,:,2)/3;
s=(I(:,:,:,1)/3).^3-(I(:,:,:,1).*I(:,:,:,2))/6+I(:,:,:,3)/2;
v=v1;

phi=acos((s./v).*sqrt((1./v)))/3;

EValue(:,:,:,1)=(I(:,:,:,1)/3+2*sqrt(v).*cos(phi));
EValue(:,:,:,2)=(I(:,:,:,1)/3-2*sqrt(v).*cos((pi/3)+phi));
EValue(:,:,:,3)=(I(:,:,:,1)/3-2*sqrt(v).*cos((pi/3)-phi));

EValue(isnan(EValue))=0;
EValue(EValue<0)=0;
clear r c v v1 s I
EValue=real(EValue);

ADC(:,:,:)=(EValue(:,:,:,1)+EValue(:,:,:,2)+EValue(:,:,:,3))/3;
Trace(:,:,:)=(EValue(:,:,:,1)+EValue(:,:,:,2)+EValue(:,:,:,3));

SqDifE=(EValue(:,:,:,1)-ADC).^2+(EValue(:,:,:,2)-ADC).^2+(EValue(:,:,:,3)-ADC).^2; 
SumSqE=EValue(:,:,:,1).^2+EValue(:,:,:,2).^2+EValue(:,:,:,3).^2;


FA=sqrt(3*(SqDifE))./sqrt(2*SumSqE);
FA(isnan(FA))=realmin;
FA(FA>1)=realmin;FA(FA<0)=realmin;
% 
ADC(find(ADC<0))=0;
Trace(find(Trace<0))=0;

'done'
