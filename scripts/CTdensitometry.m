function T = CTlung_Unreg(tag,img,voxvol,label)

T = [];

% Filter the data
if ismember(tag,{'exp','ins'})
    % Find image slices
    ind = find(std(single(img),1,[1 2])~=0)';
    if length(ind) == size(img,3)
        img = medfilt3(img);
    else
        for i = ind
            img(:,:,i) = medfilt2(img(:,:,i));
        end
    end
else
    warning('Invalid input data type: %s',tag);
    return;
end

T = lobeLoop(label,@(mask,tag,img,voxvol)unreg_sub(mask,tag,img,voxvol),...
    tag,img,voxvol);

function T = unreg_sub(mask,tag,img,voxvol)

if strcmp(tag,'exp')
    vars = {'Exp_Vol',      'double';...
            'Exp_HU',       'double';...
            'Exp_856',      'double'};
else
    vars = {'Ins_Vol',      'double';...
            'Ins_HU',       'double';...
            'Ins_950',      'double';...
            'Ins_810',      'double';...
            'Ins_810low',   'double';...
            'Ins_500',      'double';...
            'Ins_GGOI',     'double';...
            'Ins_FIBI',     'double'};
end
T = table('Size',[1,size(vars,1)],'VariableTypes',vars(:,2)','VariableNames',vars(:,1)');

nvox = nnz(mask);
maskvals = img(mask);
T{1,1} = nvox * voxvol / 1e6;
T{1,2} = mean(maskvals);
if strcmp(tag,'exp')
    T.Exp_856 = 100 * nnz(maskvals < -856) / nvox;
else % 'ins'
    T.Ins_950 = 100 * nnz(maskvals < -950) / nvox;
    T.Ins_810 = 100 * nnz((maskvals >= -810) & (maskvals < -250)) / nvox;
    T.Ins_810low = 100 * nnz((maskvals >= -810) & (maskvals < -500)) / nvox;
    T.Ins_500 = 100 * nnz((maskvals >= -500) & (maskvals < -0)) / nvox;
    T.Ins_GGOI = 100 * nnz((maskvals >= -810) & (maskvals < -700)) / nvox;
    T.Ins_FIBI = 100 * nnz((maskvals >= -700) & (maskvals < -250)) / nvox;
end
