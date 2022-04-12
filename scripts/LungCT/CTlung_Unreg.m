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

lobe = getLobeTags(label);
n = numel(lobe);
% ulab = unique(label(label>0));
% n = numel(ulab);
for i = 0:n
    if i==0
        % Use whole lung mask
        tmask = logical(label);
        S(i+1).tag = 'WholeLung';
    elseif n<1
        break;
    else
        % Use sub-region in label
        tmask = label==lobe(i).val;
        S(i+1).tag = lobe(i).name;
    end
    nvox = nnz(tmask);
    maskvals = img(tmask);
    S(i+1).vol = nvox*voxvol/1e6;
    S(i+1).mean = mean(maskvals);
    if strcmp(tag,'exp')
        S(i+1).exp856 = 100 * nnz(maskvals < -856) / nvox;
        if nargin==5 && ~isempty(atMap)
            S(i+1).SNpct = 100 * nnz(atMap(tmask)) / nvox;
            S(i+1).SNmean = mean(img(logical(atMap)));
        end
    else % 'ins'
        S(i+1).ins950 = 100 * nnz(maskvals < -950) / nvox;
        S(i+1).ins810 = 100 * nnz((maskvals >= -810) & (maskvals < -250)) / nvox;
        S(i+1).ins810low = 100 * nnz((maskvals >= -810) & (maskvals < -500)) / nvox;
        S(i+1).ins500 = 100 * nnz((maskvals >= -500) & (maskvals < -0)) / nvox;
    end
end