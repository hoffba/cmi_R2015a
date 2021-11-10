function j = Main_Pipeline_Data4DLsegDev(fpath)

%% Find files
fnames = dir(fullfile(fpath,'*.nii.gz'));
fnames = {fnames.name};
fnbase = extractBefore(fnames,'.');
id = extractBefore(fnbase,'_');
uid = unique(id);

%% Loop over Subject ID
nid = numel(uid);
for i = 1:nid
    procdir = fullfile(fpath,id{i});
    
    %% Move and rename files into SubjectID folders
    if ~isfolder(procdir)
        mkdir(procdir);
    end
                    
    [fname,ext] = findFile(fnames,{uid{i},'EXP','.ct'});
    if isempty(fname)
        fprintf('%s - missing EXP image.\n',uid{i});
    else
        fn_exp = fullfile(procdir,[uid{i},'.exp',ext]);
        movefile(fullfile(fpath,fname),fn_exp);
        
        [fname,ext] = findFile(fnames,{uid{i},'INSP','.ct'});
        if isempty(fname)
            fprintf('%s - missing EXP image.\n',uid{i});
        else
            fn_ins = fullfile(procdir,[uid{i},'.ins',ext]);
            movefile(fullfile(fpath,fname),fn_ins);
            
            [fname,ext] = findFile(fnames,{uid{i},'EXP','.lobe_segmentation'});
            if ~isempty(fname)
                fn_explabel = fullfile(procdir,[uid{i},'.exp.label',ext]);
                movefile(fullfile(fpath,fname),fn_explabel);
            end
                
            [fname,ext] = findFile(fnames,{uid{i},'INSP','.lobe_segmentation'});
            if ~isempty(fname)
                fn_inslabel = fullfile(procdir,[uid{i},'.ins.label',ext]);
                movefile(fullfile(fpath,fname),fn_inslabel);
            end
            
            % Start pipeline batch job
            fprintf('%u - starting Main_Pipeline_sub as batch job: #',uid{i});
            j(i) = batch(@Main_Pipeline_sub,1,{fn_exp,fn_ins,procdir});
            fprintf('%u\n',j(i).ID);
            
        end
    end
end

function [fname,ext] = findFile(fnames,filt)
ind = true(numel(filt),1);
for i = 1:numel(filt)
    ind = ind & contains(fnames,filt{i});
end 
fname = fnames{ind};
[~,tname,ext] = fileparts(fname);
if strcmp(ext,'.gz')
    [~,~,ext2] = fileparts(tname);
    ext = [ext2,ext];
end
