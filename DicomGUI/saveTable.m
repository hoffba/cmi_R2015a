function saveTable( h, ~, ~)
    g = guihandles(h);
    table = get(g.dcmtable, 'Data');
    f = cell2mat(table(:,1));
    num = size(f,1);
    d = uigetdir(pwd, 'Choose which directory to store files');
    pdir = '';
    %Parallel processed version, needs optimization toolbox
    %parfor x = 1:num
    %    if f(x,1)
    %    end
    %end
    for x = 1:num
        if f(x,1)
            dname = table{x,2};
            if ~strcmp(pdir, dname)
                list = dir2cell(dname);
                list = list(3:end, :);
                for i=1:length(list)
                    list(i,:) = {[dname '/' list{i,:}]};
                end
                %isolate dicoms before loading
                temp = list;
                m = strfind(temp, '.1');
                m = cellfun(@isempty, m);
                temp(m) = [];
                if isempty(temp)
                    temp = list;
                    m = strfind(temp, '.dcm');
                    m = cellfun(@isempty, m);
                    temp(m) = [];
                    if isempty(temp)
                        temp = list;
                        m = strfind(temp, '.');
                        m = cellfun(@isempty, m);
                        m = ~m;
                        temp(m) = [];
                    end
                end
                if ~isempty(temp)
                    list = temp;
                    [img, label, fov, info] = cmi_load(dname, [], list);
                    if ~exist(fullfile(d, 'Processed'),'dir');
                        mkdir(fullfile(d, 'Processed'));
                    end
                    %Isolate variables for different CSV formats
                    s = table(1,:);
                    s{1,1} = 'a';
                    nameind = strfind(s, 'SeriesDescription');
                    nameind = ~cellfun(@isempty, nameind);
                    nameind = find(nameind == 1);
                    numind = strfind(s, 'SeriesNumber');
                    numind = ~cellfun(@isempty, numind);
                    numind = find(numind == 1);
                    patind = strfind(s, 'PatientID');
                    patind = ~cellfun(@isempty, patind);
                    patind = find(patind == 1);
                
                    m = table{x,nameind};
                    m = strrep(m, ' ', '_');
                    m = strrep(m, '/', '_');
                    if ~isempty(strfind(d, 'Processed'))
                        fname = fullfile(d, [num2str(table{x,patind}) '_' m '.mhd']);
                    else
                        fname = fullfile(d, 'Processed', [num2str(table{x,patind}) '_' m '_' num2str(table{x,numind}) '.mhd']);    
                    end
                    saveMHD(fname, img, label, fov, info);
                else
                    disp(['Folder ' num2str(x) ' contained no readable Dicoms']) 
                end
            end
            pdir = dname;
        end
    end
    
    disp('Conversion Completed');
end