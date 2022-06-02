function T = tabulatePRM(mask,prm,vals,tag)
%     if flag % 10-color
%         vals = 1:10;
%         tag = cellfun(@num2str,num2cell(vals),'UniformOutput',false);
%     else % 5-color
%         vals = 1:5;
%         tag = {'Norm', 'fSAD', 'Emph', 'PD', 'NS'};
%     end
    nv = numel(vals);
    vname = strcat('PRM_',tag);
    T = table('Size',[1,nv],'VariableTypes',repmat({'double'},1,nv),...
        'VariableNames',vname);
    np = nnz(mask);
    for i = 1:numel(vals)
        T.(vname{i}) = nnz(prm(mask)==i)/np*100;
    end