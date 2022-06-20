function T = Collate_Pipeline_Results(fn_results,save_path,varnames)
% Collates pipeline results from .csv files in fn_results
% (optional) populate collated results using only variables listed in hdr

T = [];

if nargin<3
    varnames = {};
end

N = numel(fn_results);
for i = 1:N
    tT = readtable(fn_results{i});
    addResultToTable(T,tT,varnames);
end


nowstring = char(datetime('now','Format','yyyyMMddHHmmss'));
writetable(T,fullfile(save_path,['Collated_Pipeline_Results_',nowstring,'.csv']))