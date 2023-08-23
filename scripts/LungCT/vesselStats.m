% Tabulate statistics for each lobe
function T = vesselStats(mask,vessels,csa)

vars = {'VOLUME_L',               'double';...
        'VESSEL_VOLUME_L',        'double';...
        'NUM_VESSELS',          'uint32';...
        'NUM_COMPONENTS',       'uint32';...
        'NUM_ENDPOINTS',        'uint32';...
        'CSA_EXP_A',            'double';...
        'CSA_EXP_B',            'double';...
        'VESSEL_VOLUME_5DOWN_L',  'uint32';...
        'VESSEL_VOLUME_5UP_L',    'uint32'};
T = table('Size',[1,size(vars,1)],'VariableTypes',vars(:,2)','VariableNames',vars(:,1)');

V = logical(vessels) .* mask; % vessels in this lobe

C = csa .* mask; % CSA of vessels in this lobe

np = nnz(mask);
voxvol = .00625^3; % liters
T.VOLUME_L = np * voxvol; % Volume in microliters, assume voxel size of (0.625mm)^3
T.VESSEL_VOLUME_L = nnz(V) * voxvol;

if any(C,'all')
    [T.NUM_VESSELS, T.NUM_COMPONENTS, T.NUM_ENDPOINTS] = CSA_size_metrics(C);
    [T.CSA_EXP_A, T.CSA_EXP_B] = CSA_metrics(C);
    N = nnz(getCSAvessVol(V,C,0,5));
    T.VESSEL_VOLUME_5DOWN_L = N * voxvol;
    T.VESSEL_VOLUME_5UP_L = (nnz(V) - N) * voxvol;
end
