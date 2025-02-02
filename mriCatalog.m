function C = mriCatalog(studydir)
% Catalog Bruker MRI study data and convert to MHD

hdr = {'Convert','Series','Protocol','Method','TR','TE','NE','Matrix','Thk','Nslc','MT','MToff'};
nlab = length(hdr);

dnames = dir(studydir);
dnames(1:2) = [];
nd = length(dnames);
C = cell(nd,nlab);
stat = false(nd,1);

% Find scan parameters
for i = 1:nd
    if dnames(i).isdir
        pnames = {fullfile(studydir,dnames(i).name,'method')};
        if exist(pnames{1},'file')
            if exist(fullfile(studydir,dnames(i).name,'acqp'),'file')
                pnames{2} = fullfile(studydir,dnames(i).name,'acqp');
            end
%             protname = '';
            p = readBrukerMRIpar(pnames);
%             if isfield(p,'ACQ_scan_name')
%                 protname = p.ACQ_scan_name;
%             end
            C(i,:) = [{true,dnames(i).name},parsePars(p,{'ACQ_scan_name','Method',...
                'PVM_RepetitionTime','PVM_EchoTime','PVM_NEchoImages','PVM_Matrix',...
                'PVM_SliceThick','PVM_SPackArrNSlices','PVM_MagTransOnOff','PVM_MagTransOffset'})];
%             C(i,:) = {true,dnames(i).name,protname,p.Method,p.PVM_RepetitionTime,p.PVM_EchoTime,...
%                 p.PVM_NEchoImages,p.PVM_Matrix,p.PVM_SliceThick,p.PVM_SPackArrNSlices,...
%                 p.PVM_MagTransOnOff,p.PVM_MagTransOffset};
            stat(i) = true;
        end
    end
end
C = C(stat,:);
[~,ord] = sort(str2double(C(:,2)));
C = C(ord,:);
saveCell2Txt([hdr;C],fullfile(studydir,'Catalog.tsv'));

% UI for converting to MHD:
ind = cellfun(@(x)~ischar(x)&&(numel(x)>1),C);
C(ind) = cellfun(@num2str,C(ind),'UniformOutput',false);
hf = figure('Position',[900,500,800,700],'MenuBar','none','ToolBar','none',...
    'Name','Close to convert selected images to MHD ...',...
    'ToolBar','none','CloseRequestFcn','uiresume');
ht = uitable(hf,'Units','normalized','Position',[0,0,1,1],'ColumnName',hdr,...
    'ColumnWidth',{50,50,100,100,50,50,50,80,50,50,50,50},...
    'ColumnEditable',[true,false(1,nlab-1)],'Data',C);
uiwait;
stat = cell2mat(ht.Data(:,1));
delete(hf);

function v = parsePars(p,fieldnames)
nf = length(fieldnames);
v = cell(1,nf);
for i = 1:nf
    if isfield(p,fieldnames{i})
        v{i} = p.(fieldnames{i});
    end
end

% Save selected images as MHD:
% ind = find(stat);
% if ~isempty(ind)
%     odir = fullfile(studydir,'MHDs');
%     if ~exist(odir,'dir')
%         mkdir(odir);
%     end
%     ni = length(ind);
%     for i = 1:ni
%         fprintf('Processing: %u of %u\n',i,ni);
%         N = C{ind(i),2};
%         protname = regexp(C{ind(i),3},'<(\w*)','tokens');
%         protname = protname{1}{1};
%         tag = strsplit(protname,'_');
%         tag = tag{end};
%         if strcmp(tag,'Localizer')
%             %skip
%         elseif strcmp(tag,'Fat')
% %             fnout = fullfile(odir,sprintf('%s_%s_TE_%.0fus',N,protname,C{ind(i),5}*1000));
% %             [img,label,fov] = readBrukerMRI(fullfile());
%         elseif strcmp(tag,'MT')
%             fnout = fullfile(odir,sprintf('%s_%s_Offset_%.0fHz.mhd',N,protname,C{ind(i),12}));
%             [img,label,fov] = readBrukerMRI(fullfile(studydir,N,'pdata','1','2dseq'));
%             saveMHD(fnout,img,label,fov);
%         end
%     end
% end
