% Tabulate statistics for each lobe
function T = vesselStats(id,ct,seg,vessels,csa)

fprintf('Tabulating vessel results:\n');

%     lobe = getLobeTags(seg);
lobe = unique(label(label>0));

nlobes = numel(lobe);
vars = {'ID',                   'cellstr';...
    'LOBE',                 'cellstr';...
    'VOLUME',               'double';...
    'VESSEL_VOLUME',        'double';...
    'PER_EMPH',             'double';...
    'NUM_VESSELS',          'uint32';...
    'NUM_COMPONENTS',       'uint32';...
    'NUM_ENDPOINTS',        'uint32';...
    'CSA_EXP_A',            'double';...
    'CSA_EXP_B',            'double';...
    'VESSEL_VOXELS_5DOWN',  'uint32';...
    'VESSEL_VOXELS_5UP',    'uint32'};
T = table('Size',[nlobes+1,size(vars,1)],'VariableTypes',vars(:,2)','VariableNames',vars(:,1)');

for i = 0:length(lobe)
    ii = i+1;
    if i==0 % Whole-lung
        mask = ismember(seg,[lobe.val]);
        T.LOBE{ii} = 'WholeLung';
    else
        mask = seg == lobe(i).val;
        T.LOBE{ii} = lobe(i).name;
    end
    np = nnz(mask);

    T.ID{ii} = id;

    V = logical(vessels) .* mask; % vessels in this lobe

    C = csa .* mask; % CSA of vessels in this lobe

    T.VOLUME(ii) = np * 0.625^3; %Convert num voxels into volume by multiplying by voxel dim
    T.VESSEL_VOLUME(ii) = nnz(V) * 0.625^3;
    T.PER_EMPH(ii) = nnz((ct < -950) & mask) / np * 100;

    fprintf('   Calculating size metrics\n');
    [T.NUM_VESSELS(ii), T.NUM_COMPONENTS(ii), T.NUM_ENDPOINTS(ii)] = CSA_size_metrics(C);
    [T.CSA_EXP_A(ii), T.CSA_EXP_B(ii)] = CSA_metrics(C);
    T.VESSEL_VOXELS_5DOWN(ii) = nnz(getCSAvessVol(V,C,0,5));
    T.VESSEL_VOXELS_5UP(ii) = nnz(V) - T.VESSEL_VOXELS_5DOWN(ii);
    %         [T.VESSEL_VOLUME_5DOWN(ii), T.NUM_VESSELS_5DOWN(ii)] = CSA_range_metrics(C < 5 & C > 0, C);
    %         [T.VESSEL_VOLUME_5UP(ii), T.NUM_VESSELS_5UP(ii)] = CSA_range_metrics(C > 5, C);

end
end
