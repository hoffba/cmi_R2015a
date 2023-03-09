function res = addResultsToTable(res,T)
% Adds new data to existing results table
% Inputs:   res = existing results table
%           T   = new results table to add

if isempty(res)
    res = T;
else
    dT = size(T);
    dres = size(res);
    fld_T = T.Properties.VariableNames;
    fld_res = res.Properties.VariableNames;

    % add keys
    key = (1:dT(1))';
    T = addvars(T,key,'Before',1);
    key = (1:dres(1))'-dres(1);
    res = addvars(res,key,'Before',1);

    % join tables
    res = outerjoin(res,T,'mergekeys',true);
    res = removevars(res,'key');

    % rearrange variables for consistent order of columns
    ind = ismember(fld_T,fld_res);
    for i = 1:numel(fld_T)
        if ~ind(i)
            fldname = fld_T{i};
            next_match = find(ind(i+1:end),1)+i;
            if ~isempty(next_match)
                res = movevars(res,fldname,'Before',fld_T{next_match});
            end
        end
    end
end


