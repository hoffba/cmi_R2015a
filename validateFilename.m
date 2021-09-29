function fname = validateFilename(fname)

% Assume directory name is valid
[fpath,fname,ext] = fileparts(fname);
fname(regexp(fname,'[^A-Za-z0-9\._]')) = [];
fname = fullfile(fpath,[fname,ext]);


