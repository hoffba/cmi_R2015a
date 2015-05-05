function BatchRun(filt, imp)
%Run batch program, once it executes get user input about which one to
%save, then run batch program to process into mhd.
%Until I get distrib comp license
%j = batch(@dcmCatalog_NonCSV, 1, {filt});
%wait(j);
%load(j, 'C');
%possible parallel processing for higher speeds?
if nargin == 1
    imp = false;
end
if imp
    [cs, path] = uigetfile('*.csv');
    C = csvimport(fullfile(path, cs));
else
    C = dcmCatalog(filt);
end
dims = size(C);
dims(2) = dims(2)+1;
B = cell(dims);
B(:,2:dims(2)) = C;
B(:,1) = {false};
h = guihandles(open('DicomTable.fig'));
set(h.goButton, 'Callback', @saveTable);
set(h.dcmtable, 'Data', B);
set(h.allButton, 'Callback', @selAll);
set(h.fieldButton, 'Callback', @fieldCheck);

end