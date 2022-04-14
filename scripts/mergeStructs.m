function s1 = mergeStructs(s1,s2)

fldnames = fieldnames(s2);
for i = numel(fldnames)
    s1.(fldnames{i}) = s2.(fldnames{i});
end