function stat = dcestuff2FLD

[~,path] = uigetfile('dce_stuff.mat','Select DCE mat:');
if ischar(path) && exist(fullfile(path,'dce_stuff.mat'),'file')
    % Load image data
    load(fullfile(path,'dce_stuff.mat'));
    vnames = whos('ser*');
    % Load scan info
    sinfofn = dir(fullfile(path,'*scaninfo*.txt'));
    if ~isempty(sinfofn)
        sinfo = readScanInfo(fullfile(path,sinfofn(1).name));
    else
        sinfo = [];
    end
    % Determine DCE images:
    if ~isempty(vnames) && ~isempty(sinfo)
        si = find(strcmp({sinfo(:).SeriesDescript},'DYN'),1,'last');
        svar = ['ser',num2str(sinfo(si).SeriesNumber)];
        if ~isempty(si) && exist(svar,'var')
            dynimg = TLCstruct2mat(eval(svar))/100;
            %tr = sinfo(si).TR;
            fov = sinfo(si).FOV;
            d = diff(eval(['cell2mat({',svar,'(1:2).loc})'])');
            fov(3) = size(dynimg,3) * sqrt(sum(d.^2));
            label = cellstr(num2str((1:size(dynimg,4))'));
            [prem,ename] = fileparts(path(1:end-1));
            [~,sname] = fileparts(prem);
            fnout = fullfile(path,[sname,'_',ename,'_',svar,'P.fld']);
            stat = saveFLD(fnout,dynimg,label,fov);
        end
    end
end

