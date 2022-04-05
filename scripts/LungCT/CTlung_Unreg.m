function S = CTlung_Unreg(tag,img,voxvol,label,atMap)

S = struct;

% Filter the data
if ismember(tag,{'exp','ins'})
    if length(find(std(single(img),1,[1 2])~=0)') == size(img,3)
        img = medfilt3(img);
    else
        img = medfilt2(img,[3,3]);
    end
else
    warning('Invalid input data type: %s',tag);
    return;
end

ulab = unique(label(label>0));
n = numel(ulab);
for i = 1:(n+(n>1))
    if i==1
        % Use whole lung mask
        tmask = logical(label);
        S(i).tag = 'WholeLung';
    else
        % Use sub-region in label
        tmask = label==ulab(i-1);
        S(i).tag = ulab(i-1);
    end
    nvox = nnz(tmask);
    maskvals = img(tmask);
    S(i).vol = nvox*voxvol/1e6;
    S(i).mean = mean(maskvals);
    if strcmp(tag,'exp')
        S(i).exp856 = 100 * nnz(maskvals < -856) / nvox;
        if nargin == 5
            S(i).SNpct = 100 * nnz(atMap(tmask)) / nvox;
            S(i).SNmean = mean(img(logical(atMap)));
        end
    else % 'ins'
        S(i).ins950 = 100 * nnz(maskvals < -950) / nvox;
        S(i).ins810 = 100 * nnz((maskvals >= -810) & (maskvals < -250)) / nvox;
        S(i).ins810low = 100 * nnz((maskvals >= -810) & (maskvals < -500)) / nvox;
        S(i).ins500 = 100 * nnz((maskvals >= -500) & (maskvals < -0)) / nvox;
    end
end