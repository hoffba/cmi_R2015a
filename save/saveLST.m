function saveLST(t,fname)
% Save simulated __pm.lst file with timestamps for Bruker CT sorting
% Inputs: t = vector of time stamps (seconds)
%         fname = full path for desired "pm.lst" file

if nargin==0
    return;
elseif nargin<2
    [fname,fpath] = uiputfile('*.lst','Save LST file as:');
    if ischar(fname)
        fname = fullfile(fpath,fname);
    else
        return
    end
end

if ~(ischar(fname) && length(fname)>6 && strcmp(fname(end-5:end),'pm.lst'))
    error('Invalid file name ... must end with "pm.lst"');
end

fid = fopen(fname,'w');
if fid
    n = length(t);
    fprintf(fid,'[PM timemarks]\rtotal=%04u\r',n);
    for i = 1:n
        % Separate time into chunks:
            % Hours
            tt(1) = floor(t(i)/3600);
            t(i) = t(i) - tt(1)*3600;
            % Minutes
            tt(2) = floor(t(i)/60);
            t(i) = t(i) - tt(2)*60;
            % Seconds
            tt(3) = floor(t(i));
            tt(4) = round((t(i) - tt(3))*1000);
            
        fprintf(fid,'%04u=%02u:%02u:%02u:%03u\r',i-1,tt);
    end
    fclose(fid);
else
    error(['Could not open file: ' fname]);
end

