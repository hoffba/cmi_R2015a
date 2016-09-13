function fldVal = getFldVal(structIn, path2fld)
% Get field value for "path2fld" (e.g., from "[s1 path] = rfindField(structIn, 'fieldName')") 
% inside structure input "structIn"
%

[dim1 dim2]= size(path2fld);
fldVal ='NA';
if (dim1>1)
    disp(' ');
    disp(['WARNING:  MULTIPLE PATHES ARE NOT ALLOWED']);
    path2fld
    disp(['PLEASE NARROW YOUR SEARCH TO A SINGLE PATH']);
    disp(' ');
    return;
else  % single path
    path_str = strSplit2(path2fld,'.');
    path_cell = cellstr(path_str);
    j1 = 1;
    lnp = length(path_cell);
    if (isequal(path_cell{lnp},''))
        disp(['EXITING:  PATH END IS NOT A LEAF']);
        return;
    else
        if (isequal(path_cell{1},'<>') && (lnp > 1))
            j1 = 2; % start from next level
        end
        fnames = fieldnames(structIn);
        % path_cell{j1}
        if (ismember(path_cell{j1}, fnames)) % correct path-start
            s1 = structIn.(path_cell{j1});
            while (j1 < lnp)
                j1=j1+1;
                s1 = s1.(path_cell{j1});
            end
            fldVal = s1;
        else
            disp('ERROR: PATH-START NOT FOUND IN STRUCTURE');
            return; 
        end
    end
end
return;
