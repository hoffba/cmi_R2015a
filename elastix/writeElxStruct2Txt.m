% Function to save Elastix parameters to .txt file
function fname = writeElxStruct2Txt(s,fname)

if (nargin==1)
    [fname,path] = uiputfile('*.txt','Save TransformParameter.txt',...
        'InitialTransform.txt');
    if fname
        fname = fullfile(path,fname);
    else
        fname = '';
    end
elseif ~ischar(fname)
    error('Invalid input: fname');
end
if ~isempty(fname)
    fid = fopen(fname,'w');
    if fid>2
        fldnames = fieldnames(s);
        for i = 1:length(fldnames)
            tval = s.(fldnames{i});
            if isnumeric(tval)
                str = sprintf(' % .8f',tval);
            elseif ischar(tval)
                str = [' "',tval,'"'];
            elseif iscellstr(tval)
                str = sprintf(' "%s"',tval{:});
            end
            fprintf(fid,['(',fldnames{i},str,')\n']);
        end
        fclose(fid);
    end
end

