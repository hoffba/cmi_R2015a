function stat = saveCell2Txt(C,fname)

stat = false;

C = cellfun(@num2str,C,'UniformOutput',false);
d = size(C);
fstr = ['%s',repmat('\t%s',1,d(2)-1),'\n'];

fid = fopen(fname,'w');
if fid>0
    for i = 1:d(1)
        fprintf(fid,fstr,C{i,:});
    end
    fclose(fid);
    stat = true;
end