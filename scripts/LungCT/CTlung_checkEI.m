function CTlung_checkEI(ct)


%% Quick segmentation for Exp/Ins Identification:
if nnz(getRespiratoryOrgans(ct.img(1).mat)) > nnz(getRespiratoryOrgans(ct.img(2).mat))
    fprintf('Swapping images due to lung volume\n');
    ct.img = ct.img([2,1]);
    ct.dcmdir = ct.dcmdir([2,1]);
end

%% Set image labels
ct.img(1).info.label = [ct.id,'_Exp'];
ct.img(2).info.label = [ct.id,'_Ins'];
ct.img(1).info.name = [ct.id,'_',ct.tp,'.exp'];
ct.img(2).info.name = [ct.id,'_',ct.tp,'.ins'];