% Tabulate statistics for each lobe
function T = vesselStats(id,ct,seg,vessels,csa)
    
    fprintf('Tabulating vessel results:\n');
    
    lobe = getLobeTags(seg);
    nlobes = numel(lobe);
    T = table('Size',[nlobes,14],...
              'VariableTypes',{'cellstr','cellstr','double','double','double',...
                               'uint32','uint32','uint32',...
                               'double','double','double',...
                               'uint32','double','uint32'},...
              'VariableNames',{'ID','LOBE','VOLUME','VESSEL_VOLUME','PER_EMPH',...
                               'NUM_VESSELS','NUM_COMPONENTS','NUM_ENDPOINTS',...
                               'CSA_EXP_A','CSA_EXP_B','VESSEL_VOLUME_5DOWN',...
                               'NUM_VESSELS_5DOWN','VESSEL_VOLUME_5UP','NUM_VESSELS_5UP'});

    for i = 1:length(lobe)
        lobe_id = lobe(i).val;
        
        T.ID{i} = id;
        T.LOBE{i} = lobe(i).name;
        
        V = vessels .* (seg == lobe_id); % vessels in this lobe
        C = csa .* (seg == lobe_id); % CSA of vessels in this lobe
        
        T.VOLUME(i) = nnz(seg == lobe_id) * 0.625^3; %Convert num voxels into volume by multiplying by voxel dim
        T.VESSEL_VOLUME(i) = nnz(V) * 0.625^3;
        T.PER_EMPH(i) = nnz((ct < -950) & (seg == lobe_id)) / nnz(seg == lobe_id)*100;
        
        fprintf('   Calculating size metrics\n');
        [T.NUM_VESSELS(i), T.NUM_COMPONENTS(i), T.NUM_ENDPOINTS(i)] = CSA_size_metrics(C);
        [T.CSA_EXP_A(i), T.CSA_EXP_B(i)] = CSA_metrics(C);
        [T.VESSEL_VOLUME_5DOWN(i), T.NUM_VESSELS_5DOWN(i)] = CSA_range_metrics(C < 5 & C > 0, C);
        [T.VESSEL_VOLUME_5UP(i), T.NUM_VESSELS_5UP(i)] = CSA_range_metrics(C > 5, C);

    end
end
