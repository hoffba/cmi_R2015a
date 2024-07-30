% getM.m - script to obtain [x y z lobeID] matrix from lobes.mat
% adapted crudely from old code of Craig
% functional form

function y = get_M(L,E,I)

L = double(L); E = double(E); I = double(I);
lobeIDs=[1 2 3 4 5 6 11 12 13 21 22 23];
data2=cat(4,E,I);
M=[];
for j=lobeIDs
    ROI=ismember(L,j);
    exp=data2(:,:,:,1);expv=exp(ROI);
    if size(data2,4)>1
        ins=data2(:,:,:,2);insv=ins(ROI);
    else
        insv=expv;
    end
    [r,c,v] = ind2sub(size(ROI),find(ROI));
    lobe = ones(size(r))*j;
    M=cat(1,M,[expv insv r c v lobe]);
end
M=[M(:,3:5) M(:,1:2) M(:,6)];
y=M;