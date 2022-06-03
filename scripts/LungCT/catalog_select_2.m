function [selected_data,opts] = catalog_select_2(varargin)
% GUI for selection of DICOM data from DICOMcatalog information
%   For use in Ins/Exp CT analysis, looking for two matching scans per time
% Input: C = table of DICOM information
%          OR = filename of catalog .csv to read in
% Output: selected_data = structure containig DICOM info for selected data
%                      .PatientID = patient name
%                      .timepoint = structure array containing timepoint info
%                                .StudyDate = string date
%                                .Scans = structure array with scan info


selected_data = [];
opts = struct('dcmpath','',...
              'sv_path','',...
              'quickreg',false,...
              'unreg',true,...
              'airway',true,...
              'scatnet',true,...
              'vessel',true,...
              'reg',true,...
              'prm',true,...
              'tprm',true,...
              'reg_seg',true);
    
gp_valid = [];
ugroups_ic = []; % unique group indices
ngroups = 0;     % number of cases available with unique ID and date
sv_path = '';
fpath = pwd;
C = [];

%% Parse inputs:
if nargin
    % Parse input directories
    ind = find(cellfun(@(x)ischar(x)&&isfolder(x),varargin),2);
    if numel(ind)
        fpath = varargin{ind(1)};
    end
    if numel(ind)>1
        sv_path = varargin{ind(2)};
    end
    
    % Parse input options
    ind = find(cellfun(@isstruct,varargin),1);
    if ind 
        fldnames = fieldnames(opts);
        for fldi = 1:numel(fldnames)
            if isfield(varargin{ind},fldnames{fldi})
                opts.(fldnames{fldi}) = varargin{ind}.(fldnames{fldi});
            end
        end
    end
end


%% Required catalog fields:
req_fields = {'SeriesDescription',...
              'PatientName',...
              'StudyDate',...
              'SeriesNumber',...
              'ConvolutionKernel',...
              'SliceThickness',...
              'Slices'};

%% Initialize the GUI:
h = initFig;
drawnow;

%% Set catalog display:
selectCatalog(opts.dcmpath);
figure(h.fig);

%% Set path for saving results:
setSavePath(opts.sv_path);

%% End script when window closes
waitfor(h.fig);

%% Set up output structure:
if ~isempty(C)
    empty_flag = false(ngroups,1);
    tagstr = cell(size(C,1),1);
    tagstr(C.Exp) = {'Exp'};
    tagstr(C.Ins) = {'Ins'};
    C = addvars(C,tagstr,'After','Ins','NewVariableNames',{'Tag'});
    tS = table2struct(C(1,3:end)); tS = tS([]);
    selected_data = struct('UMlabel',cell(1,ngroups),'StudyDate',cell(1,ngroups),'Scans',cell(1,ngroups));
    for ig = 1:ngroups
        g_ind = ugroups_ic==ig;
        tC = C( g_ind & (C.Exp | C.Ins),:);
        if isempty(tC)
            empty_flag(ig) = true;
        else
            selected_data(ig).UMlabel = tC.UMlabel{1};
            selected_data(ig).StudyDate = tC.StudyDate{1};
            selected_data(ig).Scans = tS;
            if any(tC.Exp)
                selected_data(ig).Scans =        table2struct(tC(tC.Exp,3:end));
            else
                selected_data(ig).Scans(1).Directory = '';
            end
            if any(tC.Ins)
                selected_data(ig).Scans(2) = table2struct(tC(tC.Ins,3:end));
            else
                selected_data(ig).Scans(2).Directory = '';
            end
        end
    end
    selected_data(empty_flag) = [];
    
end

%% Set up callbacks:
    function h = initFig
        % Set up figure:
        scsz = get(0,'screensize'); % screen size
        fwidth = scsz(3)/2;
        fheight = scsz(4)*2/3;
        gap = 5;
        h.fig = uifigure('Position',round([(scsz(3)-fwidth)/2, scsz(4)/6, fwidth, fheight]),...
            'Name','Select DICOMcatalog data for processing:');%,'CloseRequestFcn',@cancel_callback);
        h.tabgroup = uitabgroup(h.fig,'Position',[0,0,fwidth,fheight]);
        h.tab1 = uitab(h.tabgroup,'Title','Scans');
        h.tab2 = uitab(h.tabgroup,'Title','Options');

        %% Tab 1: tables for selecting images, and execution buttons
        h.table_select = uitable(h.tab1,'Position',[ 1 , 20 , fwidth , fheight-130 ],'CellEditCallback',@editCell);
        h.table_filter = uitable(h.tab1,'Position',[ 1 , fheight-115 , fwidth , 90 ],'CellEditCallback',@applyFilter);

        uibutton(h.tab1,'Position',[1,1,50,20],'Text','Clear','BackgroundColor','blue','ButtonPushedFcn',@clearTags);
        uibutton(h.tab1,'Position',[fwidth-50,1,50,20],...
            'Text','Done','BackgroundColor','green','ButtonPushedFcn',@done_callback);
        uibutton(h.tab1,'Position',[fwidth-100-gap,1,50,20],...
            'Text','Cancel','BackgroundColor','red','ButtonPushedFcn',@cancel_callback);

        %% Tab 2: Settings and options
        uibutton(h.tab2,  'Position',[5,   fheight-50, 100,        20],'Text','Select Catalog:','ButtonPushedFcn',@selectCatalog);
        h.text_cat = uitextarea(h.tab2,'Position',[110, fheight-50, fwidth-115, 20],'Editable',0);
        uibutton(h.tab2,  'Position',[5,   fheight-75, 100,        20],'Text','Save To:','ButtonPushedFcn',@setSavePath);
        h.text_save = uitextarea(h.tab2,'Position',[110, fheight-75, fwidth-115, 20],'Editable',0);

        h.panel_modules = uipanel(h.tab2,'Position',[5, fheight-305, 215, 225],'Title','Modules');
        uibutton(h.panel_modules,  'Position',[5,   185, 100, 20],'Text','Select All','ButtonPushedFcn',@selectAll);
        uibutton(h.panel_modules,  'Position',[110, 185, 100, 20],'Text','Clear All','ButtonPushedFcn',@clearAll);
        h.unreg =   uicheckbox(h.panel_modules,'Position',[5,   160, 200, 20],'Text','Unreg','Value',opts.unreg);
        h.airway =  uicheckbox(h.panel_modules,'Position',[5,   135, 200, 20],'Text','Airways','Value',opts.airway);
        h.scatnet = uicheckbox(h.panel_modules,'Position',[5,   110, 200, 20],'Text','ScatNet','Value',opts.scatnet);
        h.vessel =  uicheckbox(h.panel_modules,'Position',[5,    85, 200, 20],'Text','Vessels','Value',opts.vessel);
        h.reg =     uicheckbox(h.panel_modules,'Position',[5,    60, 200, 20],'Text','Registration','Value',opts.reg);
        h.prm =     uicheckbox(h.panel_modules,'Position',[5,    35, 200, 20],'Text','PRM','Value',opts.prm);
        h.tprm =    uicheckbox(h.panel_modules,'Position',[5,    10, 200, 20],'Text','tPRM','Value',opts.tprm);

        h.panel_reg = uipanel(h.tab2,'Position',[225, fheight-305, 215, 225],'Title','Reg Options');
        h.reg_seg = uicheckbox(h.panel_reg,'Position',[5, 185, 200, 20],'Text','Tranform Segmentation','Value',opts.reg_seg);
        h.quickreg = uicheckbox(h.panel_reg,'Position',[5,160,200,20],'Text','Quick Registration','Value',opts.quickreg);

    end

    function applyFilter(~,eventdata)
        nscans = size(h.table_select.Data,1);
        
        if nargin==0 % filter all
            rowi = 1:2;
        elseif nargin==2 % filter single
            rowi = eventdata.Indices(1);
        end
        
        % Set tags:
        for i = rowi
            if any(~cellfun(@isempty,table2cell(h.table_filter.Data(i,2:end))))
                TF = true(nscans,1);
                filtnames = h.table_filter.ColumnName(2:end);
                filtnames(cellfun(@isempty,table2cell(h.table_filter.Data(i,2:end)))) = [];
                for j = 1:numel(filtnames)
                    filtval = h.table_filter.Data.(filtnames{j}){i};
                    if iscell(h.table_select.Data.(filtnames{j}))
                        % String comparison
                        TF = TF & contains(C.(filtnames{j}),filtval);
                    elseif isnumeric(h.table_select.Data.(filtnames{j}))
                        % Numeric comparison
                        TF = TF & (C.(filtnames{j})==str2double(filtval));
                    end
                end
            else
                TF = false(nscans,1);
            end
            h.table_select.Data.(h.table_filter.Data{i,1}{1}) = TF;
        end

        checkValid;
    end

    function editCell(hObject,eventdata)
        if ismember(eventdata.Indices(2),[1,2])
            % Edited EXP/INS Tag
            checkValid(hObject);
        elseif isempty(regexp(eventdata.NewData , '[/\*:?"<>|]', 'once')) % Edited UMlabel
            % Set ID for all scans in this case group:
            hObject.Data(ugroups_ic==ugroups_ic(eventdata.Indices(1)),3) = {eventdata.NewData};
        else
            warning('Invalid PatientID, try again.');
            hObject.Data{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.PreviousData;
        end
    end

    function checkValid(~)
        bgs_white = '#ffffff';
        bgs_blue =  '#bbdffb';
        bgs_green = '#99edc3';
        bgs_red =   '#ffc4c4';
        for i = 1:ngroups
            irow = find(ugroups_ic==i);
            n_exp = nnz(h.table_select.Data.Exp(irow));
            n_ins = nnz(h.table_select.Data.Ins(irow));
            gp_valid(i) = n_exp + n_ins + (n_exp>1) + (n_ins>1);
            switch gp_valid(i)
                case 0 % Nothing selected
                    tcolor = bgs_white;
                case 1 % Either Exp OR Ins selected
                    tcolor = bgs_blue;
                case 2 % Both Exp AND Ins selected
                    tcolor = bgs_green;
                otherwise % Invalid number of Exp or Ins selected for case
                    tcolor = bgs_red;
            end
            if ~mod(i,2) % Even numbered groups, make color darker
                hsv_color = rgb2hsv([ double(hex2dec(tcolor(2:3))),...
                                      double(hex2dec(tcolor(4:5))),...
                                      double(hex2dec(tcolor(6:7)))] /255);
                hsv_color(3) = hsv_color(3)*0.9;
                tcolor = hsv2rgb(hsv_color);
            end
            addStyle(h.table_select,uistyle('BackgroundColor',tcolor),'row',irow);
        end
    end
    function clearTags(~,~)
        h.table_select.Data(:,1:2) = {false};
        checkValid;
    end
    function done_callback(~,~)
        % Check that each timepoint has Ins/Exp selected (and only one of each)
        if any(gp_valid>2)
            warning('No case can have more than one Exp/Ins selected.');
        else
            C = h.table_select.Data;
            
            % Set up options structure
            opts.savepath = sv_path;
            flds = {'unreg','airway','scatnet','vessel','reg','prm','tprm','quickreg','reg_seg'};
            for iflds = 1:numel(flds)
                opts.(flds{iflds}) = h.(flds{iflds}).Value;
            end
            
            h.fig.delete;
        end
    end
    function cancel_callback(~,~)
        C = [];
        opts = [];
        h.fig.delete;
    end
    function setOpts(hObject,~)
        
    end
    function selectCatalog(str,~)
        if nargin==0 || isempty(str) || ~ischar(str)
            str = uigetdir('Select folder for processing.');
        elseif istable(str)
            C = str;
            str = '';
        end
        if isfolder(str)
            % Find existing DICOMcatalog files:
            fnames = dir(fullfile(fpath,'*DICOMcatalog*.csv'));
            answer = listdlg('PromptString','Select a catalog:',...
                'ListString',{'* Create New *',fnames.name});
            if isempty(answer)
                str = '';
            elseif answer == 1
                % Generate DICOM catalog:
                C = dcmCatalog(fpath,str);
                str = '';
            elseif answer > 1
                str = fullfile(fnames(answer-1).folder,fnames(answer-1).name);
            end
        end
        if exist(str,'file')
            iopt = detectImportOptions(str);
            C = readtable(str,iopt);
        end
        if ~isempty(C) % Set up the tables and data:
            %% Validate table input:
            colnames = C.Properties.VariableNames;
            i_member = ismember(req_fields,colnames);
            if ~all(ismember(req_fields,colnames))
                error(['Missing required fields: ',strjoin(req_fields(~i_member),', ')]);
            end
            
            %% Set up display data:
            % Re-order columns for easier reading:
            C = movevars(C,req_fields,'Before',1);

            % Add Tag columns:
            nscans = size(C,1);
            UMlabel = C.PatientName;
            ind = cellfun(@(x)isempty(x)||strcmp(x,'Anonymous'),UMlabel);
            UMlabel(ind) = C.PatientID(ind);
            C = addvars(C,false(nscans,1),false(nscans,1),UMlabel,'Before',1,'NewVariableNames',{'Exp','Ins','UMlabel'});
            colnames = C.Properties.VariableNames;
            nfields = length(colnames);

            % Remove empty data and sort array
            C(cellfun(@isempty,C.UMlabel),:) = [];
            C = sortrows(C,{'UMlabel','StudyDate','StudyID','SeriesNumber'});

            % Find groups of scans with unique PatientName and StudyDate
            if isnumeric(C.StudyDate)
                C.StudyDate = cellfun(@num2str,num2cell(C.StudyDate),'UniformOutput',false);
            end
            [~,~,ugroups_ic] = unique(strcat(C.PatientName,C.StudyDate,C.StudyID));
            ngroups = max(ugroups_ic);

            gp_valid = zeros(ngroups,1);
            
            %% Set table data:
            colEdit = false(1,nfields); colEdit([1,2,3]) = true;
            colWidth = repmat({'auto'},1,nfields); colWidth(1:2) = {40,40};
            set(h.table_select,'Data',C,'ColumnEditable',colEdit,'ColumnWidth',colWidth);
            
            %% Set filter table:
            nv = size(C,2)-3;
            fC = table('Size',[2,nv],...
                'VariableTypes',repmat({'cellstr'},1,nv),...
                'VariableNames',C.Properties.VariableNames(4:end));
            fC = addvars(fC,{'Exp';'Ins'},'Before',1,'NewVariableNames',{'Tag'});
            set(h.table_filter,'Data',fC,'ColumnEditable',[false,true(1,nv)],'ColumnWidth',[{80},colWidth(3:end)]);
            
            pause(1);
            applyFilter;
        end
    end
    function setSavePath(tpath,~)
        if nargin==0 || isempty(tpath) || ~ischar(tpath)
            tpath = uigetdir('Select folder for saving results:');
        end
        if ~isempty(tpath)
            sv_path = tpath;
            h.text_save.Value = sv_path;
        end
    end
    function selectAll(~)
        set([h.unreg,h.airway,h.scatnet,h.vessel,h.reg,h.prm,h.tprm],'Value',1);
    end
    function clearAll(~)
        set([h.unreg,h.airway,h.scatnet,h.vessel,h.reg,h.prm,h.tprm],'Value',0);
    end
end
