function T = addResultToTable(T,tT,varnames)
% Add single-row table result (tT) to table T

    if size(tT,1)>1
        warning('Result to add must be single-row table.');
        return;
    end

    if ~isempty(varnames)
        tT = tT(:,ismember(tT.Properties.VariableNames,varnames));
    end
    
    if isempty(T)
        T = tT;
        vtypes = varfun(@class,T,'OutputFormat','cell');
        defvals = cell(1,size(vtypes,2));
        defvals(strcmp(vtypes,'cell')) = {''};
        defvals(strcmp(vtypes,'double')) = {nan};
        T = addprop(T,{'defvals'},{'variable'});
        T.Properties.CustomProperties.defvals = defvals;
    else
        % Add new row of defaults
        T = [T;T.Properties.CustomProperties.defvals];
        for j = 1:size(tT,2) % Loop over varables
            vname = tT.Properties.VariableNames{j};
            if ~ismember(vname,T.Properties.VariableNames)
                % Determine variable type and default value
                vtype = class(tT.(vname));
                switch vtype
                    case 'cell'
                        defval = {''};
                    case 'double'
                        defval = nan;
                    otherwise
                        defval = {};
                end

                % Add variable to table T
                T = addvars(T,repmat(defval,size(T,1),1),...
                    'After',tT.Properties.VariableNames(j-1),...
                    'NewVariableNames',vname);
            end
            T.(vname)(end) = tT.(vname);
        end
    end
