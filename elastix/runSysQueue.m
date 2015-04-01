function runSysQueue(qfname)
% Runs system commands saved in queue file

stat = true;
while stat && exist(qfname,'file')
    fid = fopen(qfname,'r');
    if fid~=-1
        str = fgetl(fid);
        if ischar(str)
            
            % Remove line from QFile:
            tqdata = fread(fid,inf,'char');
            fclose(fid);
            fid = fopen(qfname,'w');
            fwrite(fid,tqdata,'char');
            fclose(fid);
            
            % Run system command:
            system(str);
            
        else % End of file --> delete
            fclose(fid);
            delete(qfname);
        end
    else
        stat = false;
    end
end