function idx_EIV = Step03_Identify_EI(data)
% Input: <structure> data(i).tag
%                           .mat
%                           .voi
%                           .info

%% Calculate lung VOI volume
voi_vol = zeros(1,length(data));
for i = 1:length(data)
    voi_vol(i) = nnz(data(i).voi.mat);
end
[~, idx_EIV] = sort(voi_vol);

%% Calculate mean HU of lung

% idx_Exp_voi = find(voi_vol == min(voi_vol));
% idx_Ins_voi = find(voi_vol ~= min(voi_vol));
% 
% if size(img{1},3) ~= size(img{2},3) % match img to voi using slice number
%     if size(img{1},3) == size(voi{idx_Exp_voi},3)
%         idx_Exp = 1;
%         idx_Ins = 2;
%     else
%         idx_Exp = 2;
%         idx_Ins = 1;
%     end
% else
%     % calculate mean of HU in Ins_VOI for both datasets. The lower HU
%     % should be Ins.
%     if mean(img{1}(logical(voi{idx_Ins_voi})),'all') > mean(img{2}(logical(voi{idx_Ins_voi})),'all')
%         idx_Exp = 1;
%         idx_Ins = 2;
%     else
%         idx_Exp = 2;
%         idx_Ins = 1;
%     end
% end
% 
% idx_EIV = [idx_Exp, idx_Exp_voi, idx_Ins, idx_Ins_voi];