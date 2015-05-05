function fieldCheck( h, ~, ~ )
g = guihandles(h);
field = inputdlg('Enter a keyword to check for', 'Check specific boxes based on a common field', 1);
if ~isempty(field)
    table = get(g.dcmtable, 'Data');
    dims = size(table);
    fieldt = cell(dims);
    if isnan(str2double(field))
        table2 = cellfun(@num2str, table, 'UniformOutput', false);
        table2 = cellfun(@strtrim, table2, 'UniformOutput', false);
        fieldt(:, :) = {field};
        log = cellfun(@strcmp, table2, fieldt);
    else
        fieldt(:,:) = {field};
        table2 = cellfun(@num2str, table, 'UniformOutput', false);
        log = cellfun(@strcmp, table2, fieldt);
    end
    [row, ~] = find(log);
    for i = 1:length(row)
        table(row(i), 1) = {true};
    end
    set(g.dcmtable, 'Data', table);
end
end
