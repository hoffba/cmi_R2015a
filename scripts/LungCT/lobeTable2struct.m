function S = lobeTable2struct(T,S,vind)
    varnames = T.Properties.VariableNames(vind);
    for i = 1:numel(varnames)
        for j = 1:size(T,1)
            S.(sprintf('%s_%s',varnames{i},T.LOBE{j})) = T.(varnames{i})(j);
        end
    end