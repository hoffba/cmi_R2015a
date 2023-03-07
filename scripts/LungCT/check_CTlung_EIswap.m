function TF = check_CTlung_EIswap(exp,exp_voxvol,ins,ins_voxvol)

% Quick segmentation for Exp/Ins Identification:
exp_seg = logical(getRespiratoryOrgans(medfilt2_3(exp)));
ins_seg = logical(getRespiratoryOrgans(medfilt2_3(ins)));

TF = (nnz(exp_seg)*exp_voxvol) > (nnz(ins_seg)*ins_voxvol) ...
    && (mean(exp(exp_seg)) < mean(ins(ins_seg)));