function fixMHDnames(fnames)
% Function to fix the MHD header tag for ElementDataFile path
% Assume the desired name is the current MHD name
% fixMHDnames
% fixMHDnames(dirname)
% fixMHDnames(fname)
% fixMHDnames({fname1,fname2,...})

if nargin==0
    [fnames,fpath] = uigetfile('*.mhd','Select MHD to correct:','MultiSelect','on');
    if ischar(fpath)
        if ischar(fnames)
            fnames = {fnames};
        end
        fnames = cellfun(@(x)fullfile(fpath,x),fnames,'UniformOutput',false);
    end
end

if ischar(fnames)
    if isfolder(fnames)
        % Search for MHDs in given directory:
        [D,F] = dirtree(fnames,'*.mhd');
        fnames = cellfun(@(x,y)cellfun(@(z)fullfile(x,z),y,...
            'UniformOutput',false)',D,F,'UniformOutput',false);
        fnames = [fnames{:}]';
    else
        fnames = {fnames};
    end
end

if iscellstr(fnames)
    for i = 1:length(fnames)
        [~,bname,ext] = fileparts(fnames{i});
        disp(bname);
        if ~strcmp(ext,'.mhd')
            warning(' !! File must be MHD.');
        elseif ~exist(fnames{i},'file')
            warning(' !! File does not exist.');
        else
            
            % Read existing text:
            fid = fopen(fnames{i},'r');
            str = fread(fid,'*char')';
            fclose(fid);
            
            % Find ElementDataFile tag and replace file name:
            str = regexprep(str,'ElementDataFile = (.*).raw',...
                ['ElementDataFile = ',bname,'.raw']);
            
            % Write out new file:
            fid = fopen(fnames{i},'w');
            fprintf(fid,'%c',str);
            fclose(fid);
            
            disp(' ... Done');
        end
    end
end