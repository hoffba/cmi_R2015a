function results = readYACTAairways(yactapath)

results = [];

fname = dir(fullfile(yactapath,'*_airway_results.csv'));

if isempty(fname)
    warning('No airways results found.');
else
    fname = fname(end).name;
    fprintf('Reading airway results from: %s\n',fname);
    
    % Read as text:
    fid = fopen(fullfile(yactapath,fname),'rt');
    str = fread(fid,'*char');
    fclose(fid);
    
    str = strsplit(str','\n')';
    
    % Airways Analysis
    ind = find(strcmp(str,'Airway Analysis Results'),1);
    if ~isempty(ind)
        ind = ind+1;
        go = true;
        while go
            ind = ind+1;
            tstr = str{ind};
            if contains(tstr,';')
                tstr = strsplit(tstr,';');
                varname = tstr(1);
                varval = {cellfun(@str2double,tstr(2:end))};
            elseif startsWith(tstr,'[Pi linear regression values whole lung ->')
                varname = strcat('Pi_regress_',{'a','slope','R^2'});
                varval = regexp(tstr,': ([\d\.-]+)','tokens');
                varval = [varval{:}];
            else
                go = false;
                varname = {};
            end
            for i = 1:numel(varname)
                results.(fixVarName(varname{i})) = varval{i};
            end
        end
    end
    
    % Initialize 7-generation table:
    genT = table('Size',[7,12],'VariableTypes',[{'uint8'},repmat({'double'},1,11)]);
    genstr = {'Lumen [mm²] vs. Generations'         , 'LumenArea';...
              'Wall [mm²] vs. Generations'          , 'WallArea';...
              'LumenAndWall [mm²] vs. Generations'  , 'LumenAndWallArea';...
              'WP [%] vs. Generations'              , 'WP';...
              'MedianMax [HU] vs. Generations'      , 'MedMaxHU';...
              'WT [mm] vs. Generations'             , 'WT';...
              'Total Diameter [mm] vs. Generations' , 'TotalDiameter'};
    for j = 1:size(genstr,1)
        ind = find(strcmp(str,genstr{j,1}),1);
        if ~isempty(ind)
            ind = ind+2;
            genT.Properties.VariableNames = strsplit(str{ind},';');
            for i = 1:7
                genT(i,:) = cellfun(@str2double,strsplit(str{ind+i},';'),'UniformOutput',false);
            end
            results.(genstr{j,2}) = genT;
        end
    end
    
    ind = find(strcmp(str,'#BronchiectasisIndex (path-based) vs. Generations'),1);
    if ~isempty(ind)
        label = strsplit(str{ind+2},';');
        ind = ind+3;
        tstr = str{ind};
        while contains(tstr,';')
            varval = strsplit(tstr,';');
            for i = 2:numel(varval)
                results.(fixVarName([varval{1},'_',label{i}])) = str2double(varval{i});
            end
            ind = ind+1;
            tstr = str{ind};
        end
    end
end

function vname = fixVarName(vname)

vname = regexprep(vname,'%','_pct');
vname = regexprep(vname,'+','plus');
vname = matlab.lang.makeValidName(vname);
