function changeNames(tdir,oldstr,newstr)

fnames = dir(tdir);
for i = 3:length(fnames)
    newname = regexprep(fnames(i).name,oldstr,newstr);
    if ~exist(fullfile(tdir,newname),'file')
        movefile(fullfile(tdir,fnames(i).name),fullfile(tdir,newname));
    elseif ~strcmp(fnames(i).name,newname)
        delete(fullfile(tdir,fnames(i).name));
    end
end
fixMHDnames(tdir);
