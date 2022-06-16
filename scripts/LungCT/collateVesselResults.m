function T = collateVesselResults(procdir)

fn = dir(fullfile(procdir,'**\*.ct.nii'));

T = [];
for i = 1:numel(fn)
    [~,ID] = fileparts(fn(i).folder);
    resname = fullfile(fn(i).folder,[ID,'_vesselMetrics.csv']);
    if exist(resname,'file')
        tT = readtable(resname);
        tT = lobeTable2struct(tT);
        tT = addvars(tT,{ID},'Before',1,'NewVariableNames',{'ID'});
    else
        % Add default row if analysis failed
        tT = table('Size',[1,1],'VariableTypes',{'cellstr'},'VariableNames',{'ID'});
        tT.ID = {ID};
    end
    T = addResultToTable(T,tT);
end

writetable(T,fullfile(procdir,'vesselSeg_Results.csv'));