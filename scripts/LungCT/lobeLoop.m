function res = lobeLoop(seg,fcnhandle,varargin)

lobe = getLobeTags(seg);
n = numel(lobe);
res = cell(n,1);
ROI = cell(1);

for i = 1:(n+3)    
    if i==1 % Whole-lung
        mask = ismember(seg,[lobe.val]);
        ROI{i} = 'WholeLung';
    elseif n<2
        break;
    elseif i==2
        mask = ismember(seg,[lobe(startsWith({lobe.name},'R')).val]);
        ROI{i} = 'RL';
    elseif i==3
        mask = ismember(seg,[lobe(startsWith({lobe.name},'L')).val]);
        ROI{i} = 'LL';
    elseif ~ismember(lobe(i-3).name,ROI)
        mask = ismember(seg,lobe(i-3).val);
        ROI{i} = lobe(i-3).name;
    else
        break
    end
    res{i} = feval(fcnhandle,mask,varargin{:});
end
ROI = ROI';

res = vertcat(res{:});
res = addvars(res,ROI,'Before',1);