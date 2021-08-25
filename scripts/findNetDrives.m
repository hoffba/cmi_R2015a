function netdrives = findNetDrives(str)

netdrives = [];
if ispc
    [~,r] = dos('net use');
    tokens = regexp(r,'([A-Z]+:)\s+([^\f\n\r\t\v]+)','tokens');
    netdrives = cell2struct(reshape([tokens{:}],2,[])',{'Drive','Remote'},2);
    if nargin
        netdrives(~contains({netdrives.Remote},str)) = [];
    end
else
    disp('findNetDrives is not defined yet for non-PC');
end