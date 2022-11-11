function P = pipeline_blood_density(exp,ins_reg,J,exp_seg)
% Calculates blood pseudo-density change between Exp/Ins
% REF: Towards local progression estimation of pulmonary emphysema using
%       CT. Staring M, et al. Medical Physics 41, 021905 (2014)

sig = 1;

% Expiration density
D_e = imgaussfilt3(exp,sig) + 1000;
D_e(D_e < 0) = 0;

% Inspiration density
D_i = imgaussfilt3(ins_reg,sig) + 1000;
D_i(D_i < 0) = 0;

% Blood density change
P = D_i .* J - D_e; 
P(~exp_seg) = 0;


