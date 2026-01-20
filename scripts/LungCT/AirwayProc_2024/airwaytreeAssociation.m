% Assigns quantitative values to distal nodes in the airway tree
function [B,Blabel] = airwaytreeAssociation(B,Blabel,N,voxsz,M,M_label,seg,r,func_flag)

if nargin < 8
    error('Insufficient input arguments');
end

if size(B,1) == 0 || size(N,1) == 0
    error('Empty branch or node matrices');
end

j_so = find(strcmp(Blabel,'Strahler'),1);
j_lobe = find(strcmp(Blabel,'Lobe'),1);
if isempty(j_so)
    error('Strahler column not found in Blabel');
end
if isempty(j_lobe)
    error('Lobe column not found in Blabel');
end

[nB,nBcol] = size(B);

if ndims(M) == 3
    nM = length(M_label);
    categorical_data = true;
    
    unique_vals = unique(M(:));
    unique_vals = unique_vals(~isnan(unique_vals));
    if max(unique_vals) > nM || min(unique_vals) < 1
        error('PRM categorical values must be in range 1 to %d', nM);
    end
    
elseif ndims(M) == 4
    nM = size(M,4);
    categorical_data = false;
    
    if length(M_label) ~= nM
        error('M_label length must match size(M,4)');
    end
    
else
    error('M must be 3D (categorical) or 4D (multi-map) matrix');
end

if nargin<9
    if islogical(M) || max(M(:)) <= 10
        func_flag = 2*ones(1,nM);
    else
        func_flag = ones(1,nM);
    end
end

orig_nBcol = nBcol;

edges = ( 1 : (max(B(:,j_so)))+1 ) -0.5;
soCnts = histcounts(B(:,j_so),edges);
soCnt = numel(soCnts);

B = [sortrows(B,j_so), zeros(nB,1)];
nBcol = nBcol + 1;
j_status = nBcol;

j_out = nBcol + (1:nM);
i_mean = func_flag==1;
i_pct = func_flag==2;
V_label = M_label;
V_label(i_mean) = strcat(V_label(i_mean),'_mean');
V_label(i_pct) = strcat(V_label(i_pct),'_pct');
Blabel = [Blabel, 'Status', V_label];

B = [B, nan(nB,nM)];

se = calcSe(voxsz,r);
broZ = round(N(:,2:4)./voxsz);
Nse = size(se,1);

%% FIX: Create lookup table from node ID to row index in N
% This handles non-sequential node IDs properly
node_ids = N(:,1);
max_node_id = max(node_ids);
nodeID_to_row = zeros(max_node_id, 1);  % Lookup table
for k = 1:size(N,1)
    nodeID_to_row(node_ids(k)) = k;
end

for i = 1:soCnts(1)
    node_index = B(i,3);  % This is the distal NODE ID
    
    %% FIX: Use lookup table to get row index instead of assuming ID = row
    if node_index < 1 || node_index > max_node_id
        B(i, j_status) = 1;
        continue;
    end
    
    node_row = nodeID_to_row(node_index);
    if node_row < 1 || node_row > size(broZ,1)
        B(i, j_status) = 1;
        continue;
    end
    
    %% FIX: Use node_row instead of node_index for array access
    index = se + [broZ(node_row,1), broZ(node_row,2), broZ(node_row,3)];
    num = 1;
    tmpData = NaN(Nse,nM);
    valid_voxel_count = 0;
    total_voxel_count = 0;
    
    for j = 1:Nse
        if index(j,1) > 0 && index(j,2) > 0 && index(j,3) > 0
            if index(j,1)<=size(M,1) && index(j,2)<=size(M,2) && index(j,3)<=size(M,3)
                total_voxel_count = total_voxel_count + 1;
                
                segflag = true;
                if ~isempty(seg)
                    segflag = logical(seg(index(j,1),index(j,2),index(j,3)));
                end
                
                if categorical_data
                    voxel_category = M(index(j,1),index(j,2),index(j,3));
                    if segflag
                        voxel_data = zeros(1, nM);
                        if ~isnan(voxel_category) && voxel_category >= 1 && voxel_category <= nM
                            voxel_data(voxel_category) = 1;
                            valid_voxel_count = valid_voxel_count + 1;
                        end
                        tmpData(num,:) = voxel_data;
                        num = num + 1;
                    end
                else
                    voxel_data = squeeze(M(index(j,1),index(j,2),index(j,3),:))';
                    if segflag
                        if ~any(isnan(voxel_data))
                            data_sum = sum(voxel_data);
                            if data_sum >= 80
                                tmpData(num,:) = voxel_data;
                                num = num + 1;
                                valid_voxel_count = valid_voxel_count + 1;
                            end
                        end
                    end
                end
            end
        end
    end
    
    tmpData = tmpData(1:(num-1),:);
    total_sampled_voxels = num - 1;
    
    if total_sampled_voxels == 0 || isempty(tmpData)
        B(i, j_status) = 1;
        continue;
    end
    
    category_counts = sum(tmpData, 1, 'omitnan');
    branch_percentages = (category_counts / total_sampled_voxels) * 100;
    
    if any(isnan(branch_percentages))
        B(i, j_status) = 1;
        continue;
    end
    
    total_valid_tissue = sum(branch_percentages);
    if total_valid_tissue < 80
        B(i, j_status) = 2;
        continue;
    end
    
    B(i, j_status) = 0;
    if ~isempty(tmpData) && ~all(isnan(tmpData(:)))
        if any(i_mean)
            B(i,j_out(i_mean)) = branch_percentages(i_mean) / 100;
        end
        if any(i_pct)
            B(i,j_out(i_pct)) = branch_percentages(i_pct);
        end
    end
end

Nso = soCnts(1);
for i = 2:soCnt
    ii = (Nso+1):(Nso+soCnts(i));
    B(ii,:) = sortrows(B(ii,:),3,'descend');
    Nso = Nso + soCnts(i);
end

k = 0;
max_recheck = min(nB, 50000);
reChk = [1 0 0 0 0 0 0] .*ones(max_recheck,1);

for i = (soCnts(1)+1):nB
    tmpi = find(B(:,2)==B(i,3));
    
    if isempty(tmpi)
        B(i, j_status) = 1;
        continue;
    end
    
    status_vals = B(tmpi, j_status);
    valid_daughters = tmpi(~isnan(status_vals) & status_vals == 0);
    
    if isempty(valid_daughters)
        B(i, j_status) = 1;
        continue;
    end
    
    tmpa = B(valid_daughters,j_out);
    
    if var(B(valid_daughters,j_lobe)) == 0
        B(i,j_lobe) = B(valid_daughters(1),j_lobe);
    end
    
    if ~isempty(tmpa) && any(~isnan(tmpa(:)))
        B(i, j_status) = 0;
        B(i,j_out) = mean(tmpa,1,"omitnan");
    else
        if k < size(reChk,1)
            reChk(k+1,1:3) = [B(i,2) B(i,3) i];
            reChk(k+1,4:3+length(valid_daughters)) = valid_daughters';
            k = k+1;
        end
    end
end

reChk = reChk(1:k,:);
iteration_count = 0;
max_iterations = min(50, ceil(nB/1000));

while k > 0 && iteration_count < max_iterations
    iteration_count = iteration_count + 1;
    resolved_count = 0;
    
    for i = 1:size(reChk,1)
        if reChk(i,7) == 1
            continue;
        end
        
        daughter_indices = reChk(i,4:6);
        daughter_indices = daughter_indices(daughter_indices > 0);
        status_vals = B(daughter_indices, j_status);
        valid_daughters = daughter_indices(~isnan(status_vals) & status_vals == 0);
        
        if ~isempty(valid_daughters)
            parent_idx = reChk(i,3);
            if var(B(valid_daughters,j_lobe)) == 0
                B(parent_idx,j_lobe) = B(valid_daughters(1),j_lobe);
            end
            B(parent_idx,j_out) = mean(B(valid_daughters,j_out),1,"omitnan");
            B(parent_idx, j_status) = 0;
            reChk(i,7) = 1;
            resolved_count = resolved_count + 1;
        end
    end
    
    if resolved_count == 0
        break;
    end
end

B = sortrows(B,1);

end
