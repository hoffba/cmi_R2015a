function T = compile_CFF_results

basepath = 'R:\CGalban_Lab\CT_Lung\clinical\data\UM\Vibha_Xplant\CFF_Trial_20200701\CT_Data\Trial_Data';
fn_save = fullfile(basepath,'CFF_qCT_Results.xlsx');
sitename = {'MM','UPenn','UToronto'};
tabname = 'qCT';

if isfile(fn_save)
    T = readtable(fn_save);
    vtype = T.Properties.VariableTypes;
    vtype(strcmp(vtype,'cell')) = {'cellstr'};
    strvars = {'ID','ROI'}; % Force these to be string
    strvars(~ismember(strvars,T.Properties.VariableNames)) = [];
    vtype(ismember(T.Properties.VariableNames,strvars)) = {'cellstr'};
    for i = 1:numel(strvars)
        if isnumeric(T.(strvars{i}))
            T.(strvars{i}) = arrayfun(@num2str,T.(strvars{i}),'UniformOutput',false);
        end
    end
else
    [T,vtype] = initTable();
end
nv = length(vtype);
strflag = ismember(T.Properties.VariableNames,strvars);

for i = 1:3
    batchname = dir(fullfile(basepath,sitename{i},'ProcessedData','batch20*'));
    for j = 1:numel(batchname)
        batchstr = extractAfter(batchname(j).name,'batch');
        fname = dir(fullfile(batchname(j).folder,batchname(j).name,'*','*_PipelineResults.csv'));
        for k = 1:numel(fname)
            % Read pipeline results
            res = readtable(fullfile(fname(k).folder,fname(k).name));

            % Check if line exists
            caserow = find( strcmp(T.ID,res.ID{1} & strcmp(T.Batch,batchstr)) ,1);
            
            if ~caserow
                t = table('Size',[1,nv],'VariableTypes',vtype,'VariableNames',T.Properties.VariableNames);
                irow = find(strcmp(res.ROI,'WholeLung'),1);
                for ivar = 1:nv
                    vname = T.Properties.VariableNames{ivar};
                    if iscell(T.(vname))
                        tval = {''};
                    elseif isnumeric(T.(vname))
                        tval = nan;
                    end
                    if strcmp(vname,'Batch')
                        tval = batchstr;
                    elseif ismember(vname,res.Properties.VariableNames)
                        if strflag(ivar)
                            if isnumeric(res.(vname))
                                tval = {num2str(res.(vname)(irow))};
                            else
                                tval = {res.(vname)(irow)};
                            end
                        else
                            tval = res.(vname)(irow);
                        end
                    end
                    t.(vname) = tval;
                end
                T = [T;t];
            end
        end
    end
end
T = removevars(T,'ROI');
writetable(T,fn_save,'Sheet',tabname);




function [T,vtype] = initTable()
var = {'ID',                    'cellstr';...
       'Batch',                 'cellstr';...
       'ROI',                   'cellstr';...
   	   'WallPct_3_8_1',         'single';...
       'Wall_pct_1',            'single';...
   	   'Pi10',                  'single';...
   	   'Pi15',                  'single';...
       'BEI',                   'single';...
       'BEI_gen',               'single';...
       'WT_seg',                'single';...
       'WT_subseg',             'single';...
       'VESSEL_VOLUME_L',       'single';...
       'VESSEL_VOLUME_5DOWN_L', 'single';...
       'VESSEL_VOLUME_5UP_L',   'single';...
       'Exp_Vol',               'single';...
       'Exp_HU',                'single';...
       'Exp_856',               'single';...
       'Ins_Vol',               'single';...
       'Ins_HU',                'single';...
       'Ins_950',               'single';...
       'Ins_810',               'single';...
       'Ins_810low',            'single';...
       'Ins_500',               'single';...
       'Ins_GGOI',              'single';...
       'Ins_FIBI',              'single';...
       'Jac_mean',              'single';...
       'Jac_var',               'single';...
       'dBlood_mean',           'single';...
       'dBlood_var',            'single';...
       'PRM_norm_pct',          'single';...
       'PRM_fsad_pct',          'single';...
       'PRM_emph_pct',          'single';...
       'PRM_pd_pct',            'single';...
       'PRM_ns_pct',            'single'};
vtype = var(:,2)';
T = table('Size',[0,size(var,1)],'VariableTypes',vtype,'VariableNames',var(:,1)');