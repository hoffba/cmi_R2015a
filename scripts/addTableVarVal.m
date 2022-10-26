   
function T = addTableVarVal(T,varargin)
    if nargin==2 && istable(varargin{1}) && any(strcmp(varargin{1}.Properties.VariableNames,'ROI'))
        vals = varargin{1};
        regionstr = vals.ROI;
        vals = removevars(vals,'ROI');
        varstr = vals.Properties.VariableNames;
    elseif nargin==4
        varstr = varargin{1};
        regionstr = varargin{2};
        vals = varargin{3};
        if ischar(regionstr)
            regionstr = {regionstr};
        end
        if ischar(varstr)
            varstr = {varstr};
        end
        if numel(varstr)<size(vals,2)
            varstr = strcat(varstr,'_',cellfun(@num2str,num2cell(1:size(vals,2)),'UniformOutput',false));
        end
    else
        return;
    end
    if numel(regionstr)==size(vals,1)
        % Table T must have RowNames set to ROI
        if isempty(T.Properties.RowNames)
            T.Properties.RowNames = T.ROI;
        end
        for i = 1:numel(varstr)
            if istable(vals)
                tvals = vals.(varstr{i});
            else
                tvals = vals(:,i);
            end
            if iscell(tvals)
                defval = {''};
            elseif isnumeric(tvals)
                defval = nan;
            else
                defval = {[]};
            end
            if ~ismember(varstr{i},T.Properties.VariableNames)
                T = addvars(T,repmat(defval,size(T,1),1),'NewVariableNames',varstr(i));
            end
            T.(varstr{i})(regionstr) = tvals;
        end
    end
    