function T = catalog_hierarchical(basepath)

ext = 'mhd';
    
if ~nargin
    basepath = pwd;
end

% Set up table
T = [];
vars = {'CaseNumber',             'uint16';...
        'Tag',                    'cellstr';...
        'UMlabel',                'cellstr';...
        'StudyDate',              'cellstr';...
        'FileName',               'cellstr';...
        'PatientName',            'cellstr';...
        'DataPath',               'cellstr';...
        'DataType',               'cellstr'};
nvars = size(vars,1);
Tdef = table('Size',[1,nvars],'VariableTypes',vars(:,2)','VariableNames',vars(:,1)');

% Set up tag translation
tags = {'Exp',      '_exp.';...
        'Ins',      '_ins.';...
        'ExpLabel', {'_explabel.','_exp_label.'};...
        'InsLabel', {'_inslabel.','_ins_label.'};...
        'InsReg',   '_ins_R.'};

case_i = 0;
D = dir(basepath); % Folders named with PatientID
for i = 3:numel(D)
    Dpath = fullfile(D(i).folder,D(i).name);
    if isfolder(Dpath)
        TP = dir(Dpath); % Subfolders named with timestamp
        for j = 3:numel(TP)
            Tpath = fullfile(TP(j).folder,TP(j).name);
            if isfolder(Tpath)
                case_i = case_i +1;
                
                % Find image files
                F = dir(fullfile(Tpath,['*.',ext]));
                
                % Find elastix reg folder
                ELX = dir(fullfile(Tpath,'elxreg_*'));
                if ~isempty(ELX)
                    for ix = 1:numel(ELX)
                        F = [F;dir(fullfile(ELX(ix).folder,ELX(ix).name,['*.',ext]))];
                    end
                end
                
                % Add scans to table
                for k = 1:numel(F)
                    tag = tags(cellfun(@(x)endsWith(F(k).name,strcat(x,ext),'IgnoreCase',true),tags(:,2)),1);
                    if isempty(tag) || startsWith(F(k).name,'.')
                        tag = {''};
                    end
                    
                    Tnew = Tdef;
                    Tnew.CaseNumber = case_i;
                    Tnew.Tag = tag;
                    Tnew.UMlabel = {D(i).name};
                    Tnew.FileName = {F(k).name};
                    Tnew.PatientName = {D(i).name};
                    Tnew.StudyDate = {TP(j).name};
                    Tnew.DataPath = {fullfile(F(k).folder,F(k).name)};
                    Tnew.DataType = {ext};
                    
                    T = [T;Tnew];
                end
            end
        end
    end
end

if ~isempty(T)
    % Save results to CSV file:
    fn_res = fullfile(procdir,sprintf('%s_SourceData.csv',ID));
    writetable(T,fn_res);
end
