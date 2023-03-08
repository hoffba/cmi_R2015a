% SCAN: CT Files and Lobe segmentation
candidates=dir('R:\CGalban_Lab\LabMembers\Bingzhao\Data\Beaumont_I2E_aligned\*warped.nii.gz'); 

% SCAN: PRM 
gt_dir = 'R:\CGalban_Lab\LabMembers\Bingzhao\Data\Beaumont\Beaumont_I2E_PRM';

% SAVE: Save the output masks
save_dir = 'R:\CGalban_Lab\LabMembers\Bingzhao\Beaumont_I2E_Save_Demo';


% evaluation variable
Eval = [];
PatiendID = string();

for i=1:size(candidates, 1)
    
    % load the three files
    imgFlNm = fullfile(candidates(i).folder, [candidates(i).name(1: end-13) 'ct.nii.gz']);
    msFlNm = fullfile(candidates(i).folder, [candidates(i).name(1: end-13) 'njh.lobe_segmentation.nii.gz']);
    imgFlNm_warp = fullfile(candidates(i).folder, candidates(i).name);
    
    if ~exist(msFlNm)
        msFlNm = fullfile(candidates(i).folder, [candidates(i).name(1: end-13) 'imbio.lobe_segmentation.nii.gz']);
    end
    
    % segmentation
    try
        [atMap] = ScatterNetAT(imgFlNm, msFlNm, 0);
        [emphMap] = ScatterNetEMPH(imgFlNm_warp, msFlNm, 0);
        fsadMap = (atMap - emphMap) > 0;
        other_at_map = (atMap - emphMap) < 0;
    catch ME
        disp([num2str(i) '/' num2str(size(candidates, 1)) '  ' imgFlNm  '  Broken, Skipped']);
        continue;
    end
    
    % fsadMap = 2, emphMap = 3, then save as one map, filename as:
    % 10017X_EXP_STD_BWH_COPD_ScatterNet_PRM.nii.gz
    saveMap = zeros(size(fsadMap));
    saveMap(emphMap > 0) = 3;
    saveMap(other_at_map > 0) = 4;
    saveMap(fsadMap > 0) = 2;
    
    % organize the file save path
    [filepath,name,ext] = fileparts(imgFlNm);
    save_name = [name(1:end-7) '_ScatterNet_PRM'];
    save_path = fullfile(save_dir, save_name);
    niftiwrite(saveMap, save_path,'Compressed',true);
    
    % Dice evaluation
    gtMask = niftiread(fullfile(gt_dir, [candidates(i).name(1: end-13) 'prm.nii.gz']));
    
    % Lobe mask [todo]
    fsadGtMask = gtMask == 3;
    eval_fsad_dice = dice(fsadGtMask, fsadMap);
    
    emphGtMask = (gtMask == 4) + (gtMask == 5);
    emphGtMask = emphGtMask > 0;
    eval_emph_dice = dice(emphGtMask, emphMap);
    
    % Compute percentage of FSAD/EMPH in the mask
    lung_mask = gtMask > 0;
    percent_fsad = sum(fsadMap(:)) / sum(lung_mask(:));
    percent_emph = sum(emphMap(:)) / sum(lung_mask(:));
    
    % [Discuss TODO]
    percent_at = (sum(fsadMap(:)) + sum(emphMap(:)) +  sum(other_at_map(:))) / sum(lung_mask(:));
    
    % Compute percentage of FSAD/EMPH in the PRM mask
    percent_fsad_prm = sum(fsadGtMask(:)) / sum(lung_mask(:));
    percent_emph_prm = sum(emphGtMask(:)) / sum(lung_mask(:));
    percent_at_prm = percent_fsad_prm + percent_emph_prm;
    
    % Visualization: off for faster speed
    disp([num2str(i) '/' num2str(size(candidates, 1)) '  ' save_name(1:6) ' --> ' 'Dice_FSAD=' num2str(eval_fsad_dice) '   | Dice_EMPH=' num2str(eval_emph_dice) ' | FSAD_PERCENT=' num2str(percent_fsad) ' | EMPH_PERCENT=' num2str(percent_emph) ' | PRM_FSAD_PERCENT=' num2str(percent_fsad_prm) ' | PRM_EMPH_PERCENT=' num2str(percent_emph_prm)]);
    
    
    % collect evaluation data
    Eval(i, 1:8) = [eval_fsad_dice, eval_emph_dice, percent_fsad, percent_fsad_prm, percent_emph, percent_emph_prm, percent_at, percent_at_prm];
    PatiendID(i) = save_name(1:6);

end


% For evaluation use
T = array2table(Eval);
T.PatiendID = PatiendID(:);
T.Properties.VariableNames(1:9) = {'fsad_dice','emph_dice', 'percent_fsad', 'percent_fsad_prm', 'percent_emph', 'percent_emph_prm', 'percent_at', 'percent_at_prm', 'sid'};
writetable(T,fullfile(save_dir, 'evaluation_bundle.csv'))


