function results = CTlung_Pipeline()


[cases,opts] = catalog_select_3;

N = numel(cases);
results = [];
if N
    % Initialize the pipeline object
    pobj = CTlungClass(opts);
    
    % Loop over cases
    for icase = 1
        pobj.setcase([cases(icase).UMlabel,'_',cases(icase).StudyDate],procdir,fn);
        res = pobj.run;
        
        % Add results to table
        results = [results;res];
    end
end