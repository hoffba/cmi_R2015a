% Tabulate statistics for each lobe
function T = vesselStats(mask,ct,vessels,csa)

vars = {'VOLUME',               'double';...
        'VESSEL_VOLUME',        'double';...
        'PER_EMPH',             'double';...
        'NUM_VESSELS',          'uint32';...
        'NUM_COMPONENTS',       'uint32';...
        'NUM_ENDPOINTS',        'uint32';...
        'CSA_EXP_A',            'double';...
        'CSA_EXP_B',            'double';...
        'VESSEL_VOXELS_5DOWN',  'uint32';...
        'VESSEL_VOXELS_5UP',    'uint32'};
T = table('Size',[1,size(vars,1)],'VariableTypes',vars(:,2)','VariableNames',vars(:,1)');

V = logical(vessels) .* mask; % vessels in this lobe

C = csa .* mask; % CSA of vessels in this lobe

np = nnz(mask);
T.VOLUME = np * 0.625^3; %Convert num voxels into volume by multiplying by voxel dim
T.VESSEL_VOLUME = nnz(V) * 0.625^3;
T.PER_EMPH = nnz((ct < -950) & mask) / np * 100;

if any(C,'all')
    [T.NUM_VESSELS, T.NUM_COMPONENTS, T.NUM_ENDPOINTS] = CSA_size_metrics(C);
    [T.CSA_EXP_A, T.CSA_EXP_B] = CSA_metrics(C);
    T.VESSEL_VOXELS_5DOWN = nnz(getCSAvessVol(V,C,0,5));
    T.VESSEL_VOXELS_5UP = nnz(V) - T.VESSEL_VOXELS_5DOWN;
%         [T.VESSEL_VOLUME_5DOWN(ii), T.NUM_VESSELS_5DOWN(ii)] = CSA_range_metrics(C < 5 & C > 0, C);
%         [T.VESSEL_VOLUME_5UP(ii), T.NUM_VESSELS_5UP(ii)] = CSA_range_metrics(C > 5, C);
end

