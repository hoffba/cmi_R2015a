function [prm5,label] = prm_convert10to5(prm10)
% function to convert 10-color Lung CT PRM to the standard 5-color

label  = {'Norm', 'fSAD', 'Emph', 'PD',       'NS'};
group = {[1,2],  3,      [4,5],  [8,9,10],   6   };
prm5 = int8(zeros(size(prm10)));
for i = 1:5
    prm5(ismember(prm10,group{i})) = i;
end