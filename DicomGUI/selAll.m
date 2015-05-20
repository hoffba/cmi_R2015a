function selAll( h, ~, ~ )
%Select all checkmarks, very vast function
g = guihandles(h);
m = get(g.dcmtable, 'Data');
m(2:end, 1) = {true};
set(g.dcmtable, 'Data', m);
end

