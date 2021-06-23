function Step04_Rename_EI(procdir,fname_Exp,fname_Exp_Label,fname_Ins,fname_Ins_Label,data)

if ~isfolder(procdir)
    mkdir(procdir);
end

if ~exist(fullfile(procdir,[fname_Exp,'.nii.gz']),'file')
    fprintf('\n   Saving file: %s\n',[fname_Exp,'.nii.gz'])
    niftiwrite(int16(data(1).img.mat),fullfile(procdir,fname_Exp),data(1).img.info,'Compressed',true);
end
if ~exist(fullfile(procdir,[fname_Ins,'.nii.gz']),'file')
    fprintf('   Saving file: %s\n',[fname_Ins,'.nii.gz'])
    niftiwrite(int16(data(2).img.mat),fname_Ins,data(2).img.info,'Compressed',true);
end
if ~exist(fullfile(procdir,[fname_Exp_Label,'.nii.gz']),'file')
    fprintf('   Saving file: %s\n',[fname_Exp_Label,'.nii.gz'])
    niftiwrite(int8(data(1).voi.mat),fname_Exp_Label,data(1).voi.info,'Compressed',true);
end
if ~exist(fullfile(procdir,[fname_Ins_Label,'.nii.gz']),'file')
    fprintf('   Saving file: %s\n',[fname_Ins_Label,'.nii.gz'])
    niftiwrite(int8(data(2).voi.mat),fname_Ins_Label,data(2).voi.info,'Compressed',true);
end

% Save pathname/filenames
% Exp_dFile = fullfile(str_home,[ID,'_Exp']);
% Ins_dFile = fullfile(str_home,[ID,'_Ins']);
% 
% if ~exist([Exp_dFile,'.nii.gz']) && ~exist([Ins_dFile,'.nii.gz'])
%     niftiwrite(int16(img{idx_Exp}),Exp_dFile,info{idx_Exp},'Compressed',true);
%     niftiwrite(int16(img{idx_Ins}),Ins_dFile,info{idx_Ins},'Compressed',true);
% end
% 
% Exp_mFile = fullfile(str_home,[ID,'_Exp_Label']);
% Ins_mFile = fullfile(str_home,[ID,'_Ins_Label']);
% 
% if ~exist([Exp_mFile,'.nii.gz'])&&~exist([Ins_mFile,'.nii.gz'])
%     if strcmp(info_voi{1},'Not Nifti')
%         info_voi{idx_Exp_voi} = info{idx_Exp};
%         info_voi{idx_Ins_voi} = info{idx_Ins};
%         
%         info_voi{idx_Exp_voi}.Datatype = 'int8'; info_voi{idx_Ins_voi}.BitsPerPixel = 8;
%         info_voi{idx_Ins_voi}.Datatype = 'int8'; info_voi{idx_Ins_voi}.BitsPerPixel = 8;
%         
%         info_voi{idx_Exp_voi}.Filename = masks{idx_Exp_voi};
%         info_voi{idx_Ins_voi}.Filename = masks{idx_Ins_voi};
%     end
%     
%     niftiwrite(int8(voi{idx_Exp_voi}),Exp_mFile,info_voi{idx_Exp_voi},'Compressed',true);
%     niftiwrite(int8(voi{idx_Ins_voi}),Ins_mFile,info_voi{idx_Ins_voi},'Compressed',true);
% end