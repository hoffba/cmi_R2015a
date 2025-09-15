function T = compile_pipeline_results(fn,varnames)

T = [];

try

    string_vars = {'ID','Exp_Source','Ins_Source','ROI'};

    % Initialize table with requested variables
    uvars = {};
    if nargin>1
        uvars = unique([string_vars,varnames],'stable');
        T = table('Size',[0,numel(uvars)],...
            'VariableTypes',[repmat({'cellstr'},1,4),repmat({'double'},1,numel(uvars)-4)],...
            'VariableNames',uvars);
    end

    % fn = dir(fullfile(basedir,'**','*_Results.csv'));
    nf = numel(fn);
    for i = 1:nf

        fname = fullfile(fn(i).folder,fn(i).name);

        fprintf('(%d of %d) : %s\n',i,nf,fname)

        opts = detectImportOptions(fname);
        t = readtable(fname,opts);

        % Force certain variables to be cellstr
        for j = 1:numel(string_vars)
            if ismember(string_vars,t.Properties.VariableNames)
                if ~iscellstr(t.(string_vars{j}))
                    t.(string_vars{j}) = repmat({''},size(t,1),1);
                end
            end
        end

        if ~istable(T)
            T = t;
        else
            T = addtotable(T,t,uvars);
        end
    end

catch err
    disp(getReport(err))
end


function T = addtotable(T,t,uvars)

    nT = size(T,1);
    nt = size(t,1);
    ind_add = nT + (1:nt);
    vT = T.Properties.VariableNames;
    vt = t.Properties.VariableNames;
    vt = vt(ismember(vt,uvars));

    % Add missing variables to T
    % ind = find(~ismember(vt,vT));
    for i = find(~ismember(vt,vT))
        if iscell(t.(vt{i}))
            tvar = repmat({''},nT,1);
        else
            tvar = nan(nT,1);
        end
        T.(vt{i}) = tvar;
    end
    
    % Add default rows to T
    vT = T.Properties.VariableNames;
    rdef = cell(1,numel(vT));
    for i = 1:numel(rdef)
        if iscell(T.(vT{i}))
            rdef{i} = '';
        else
            rdef{i} = nan;
        end
    end
    T = [T;repmat(rdef,nt,1)];

    % Add values in t to new rows in T
    for i = 1:numel(vt)
        T.(vt{i})(ind_add) = t.(vt{i});
    end












