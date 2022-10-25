function fixTransformParameter(fn)

go = true;
while go
    if ischar(fn) && endsWith(fn,'.txt') && exist(fn,'file')
        % Read TP file:
        fid  = fopen(fn,'r');
        str = fread(fid,'*char')';
        fclose(fid);
        
        % Find Initial TP filename:
        tok = regexp(str,'InitialTransformParametersFileName "([^)]+)")','tokens');
        full_tp = tok{1}{1};
        
        if strcmp(full_tp,'NoInitialTransform')
            go = false;
        else
            % Replace full path with filename
            [~,tp,ext] = fileparts(full_tp);
            str = regexprep(str,full_tp,[tp,ext]);
            fid = fopen(fn,'w');
            fprintf(fid,'%s',str);
            fclose(fid);
            
            fn = fullfile(fileparts(fn),[tp,ext]);
        end
    end
end