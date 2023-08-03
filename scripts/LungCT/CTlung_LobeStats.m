function T = CTlung_LobeStats(seg,str,stats,varargin)
% seg =     lobe segmentation map
% str =     name of stats to tabulate
% stats =   statistics to run (e.g. 'mean', 'pct', ...)
% label =   optional, additional logical/categorical map
% img =     optional, additional image

lobe = getLobeTags(seg);
n = numel(lobe);

% Validate inputs
if ischar(stats)
    stats = {stats};
end
stats(~ismember(stats,{'mean','var','pct','prm'})) = [];
label = true;
labl_vals = true;
img = [];
vname = string(str);
for i = 1:numel(varargin)
    if ~isempty(varargin{i})
        if iscategorical(varargin{i})
            label = varargin{i};
            labl_vals = categories(label);
            vname = vname + '_' + labl_vals;
        elseif islogical(varargin{i})
            label = varargin{i};
        elseif isnumeric(varargin{i})
            img = varargin{i};
        end
    end
end
vname = vname + '_' + stats;

% Determine regions to analyze
ROI = {};
if n
    ROI = {'WholeLung', [lobe.val]};
end
if n>1
    ROI = [ROI;...
           {'RL', [lobe(startsWith({lobe.name},'R')).val]};...
           {'LL', [lobe(startsWith({lobe.name},'L')).val]}];
end
lobe = lobe(~ismember({lobe.name},ROI(:,1)));
if ~isempty(lobe)
    ROI = [ROI;...
           [{lobe.name}' , {lobe.val}']];
end

% if startsWith(str,'prm')
%     if strcmp(str,'prm10') % 10-color
%         labl_vals = 1:10;
%         tag = string(labl_vals);
%     else % 5-color
%         labl_vals = 1:5;
%         tag = ["Norm", "fSAD", "Emph", "PD", "NS"];
%     end
%     vname = "PRM_" + tag' + '_' + stats;
% else
%     labl_vals = 1;
%     vname = string(str) + "_" + stats;
% end

% Initialize results table
nr = size(ROI,1);
nv = numel(vname);
T = table('Size',[nr,nv],'VariableTypes',repmat({'double'},1,nv),'VariableNames',vname(:));

% Loop over lung lobes (including R/L and Whole-Lung)
for irow = 1:nr
    roimask = ismember(seg,ROI{irow,2});
    np = nnz(roimask);
    % Loop over label values
    for ilabl = 1:numel(labl_vals)
        if islogical(label)
            mask = roimask & label;
        else
            mask = roimask & label==labl_vals{ilabl};
        end
        % Generate table of stats for this mask
        for istat = 1:numel(stats)
            switch stats{istat}
                case 'mean'
                    if ~isempty(img)
                        T.(vname{ilabl,istat})(irow) = mean(img(mask));
                    end
                case 'var'
                    if ~isempty(img)
                        T.(vname{ilabl,istat})(irow) = var(img(mask));
                    end
                case 'pct'
                    T.(vname{ilabl,istat})(irow) = nnz(mask)/np*100;
            end
        end
    end
end
T = addvars(T,ROI(:,1),'Before',1,'NewVariableNames',{'ROI'});

