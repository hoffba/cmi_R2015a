function res = lobeLoop(seg,fcnhandle,varargin)

lobe = getLobeTags(seg);
n = numel(lobe);
res = cell(n,1);
LOBE = cell(n,1);

for i = 1:(n+3)    
    if i==1 % Whole-lung
        mask = ismember(seg,[lobe.val]);
        LOBE{i} = 'WholeLung';
    elseif n<2
        break;
    elseif i==2
        mask = ismember(seg,[lobe(startsWith({lobe.name},'R')).val]);
        LOBE{i} = 'RL';
    elseif i==3
        mask = ismember(seg,[lobe(startsWith({lobe.name},'L')).val]);
        LOBE{i} = 'LL';
    elseif ~ismember(lobe(i-3).name,LOBE)
        mask = seg == lobe(i-3).val;
        LOBE{i} = lobe(i-3).name;
    else
        break
    end
    res{i} = feval(fcnhandle,mask,varargin{:});
end

res = vertcat(res{:});
res = addvars(res,LOBE,'Before',1);