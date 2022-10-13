function [selected_data,opts] = catalog_select_3(varargin)
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
opts = struct('par_size',5,...
              'username','',...
              'dcm_path','',...
              'save_path','',...
              'quickreg',false,...
              'unreg',true,...
              'airway',true,...
              'scatnet',true,...
              'vessel',true,...
              'reg',true,...
              'prm',true,...
              'tprm',true,...
              'reg_seg',false,...
              'orient_check',true,...
              'jac',true,...
              'jacmat',false,...
              'def',false,...
              'Tvoi',true,...
              'Tsurf',false);
          
C = [];          % full catalog table

% Scan groups
ugroups_ic = []; % unique group indices
ngroups = 0;     % number of cases available with unique ID and date
gp_valid = [];   

% Page groups
gp_per_page = 5; % Number of groups to show per page
pageN = 0;       % number of pages in catalog
curr_page = 0;   % current page number
page_ic = [];    % current catalog indices for displaying (only shows 5 cases at once, per page, to limit lag)

% Set colors:
bgs_white = '#ffffff';
bgs_blue =  '#bbdffb';
bgs_green = '#99edc3';
bgs_red =   '#ffc4c4';

%% Parse inputs:
validateInputs(varargin);

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

if ~isempty(opts.dcm_path)
    selectCatalog(opts.dcm_path);
end
if ~isempty(opts.save_path)
    setSavePath(opts.save_path);
end

%% End script when window closes
waitfor(h.fig);

%% Gather outputs
if ~isempty(C)
    %% Determine tags:
    tagstr = repmat({''},size(C,1),1);
    tagstr(C.Exp) = {'Exp'};
    tagstr(C.Ins) = {'Ins'};
    C = addvars(C,ugroups_ic,tagstr,'Before',1,'NewVariableNames',{'CaseNumber','Tag'});
    
    if isa(C.UMlabel,'double') % CJG changed b/c error for double in lin 78
        C.UMlabel = cellstr(strcat('ID',int2str(C.UMlabel)));
    end

    %% Set up output structure:
    empty_flag = false(ngroups,1);
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
    
    %% Save selections to catalog file:
    fname = opts.dcm_path;
    if ~isempty(fname)
        C = removevars(C,{'Ins','Exp'});
        if isfolder(fname)
            fname = fullfile(fname,'DICOMcatalog_select.csv');
        end
        try
            writetable(C,fname);
        catch
            fname = fullfile(fileparts(fname),sprintf('DICOMcatalog_select_%s.csv',datestr(datetime('now'),'yyyymmddHHMMSS')));
            writetable(C,fname);
        end
    end

end

%% Set up callbacks:
    function validateInputs(inputs)
        p = inputParser;
        addParameter(p,'dcm_path','',@(x)ischar(x)&&(isempty(x)||isfolder(x)));
        addParameter(p,'save_path','',@(x)ischar(x)&&(isempty(x)||isfolder(x)));
        addParameter(p,'opts',opts,@isstruct);
        parse(p,inputs{:});
        opts.dcm_path = p.Results.dcm_path;
        opts.save_path = p.Results.save_path;
        flds = fieldnames(p.Results.opts);
        for i = 1:numel(flds)
            if isfield(opts,flds{i})
                opts.(flds{i}) = p.Results.opts.(flds{i});
            end
        end
    end
    function h = initFig
        % Set up figure:
        scsz = get(0,'screensize'); % screen size
        fwidth = scsz(3)/2;
        fheight = scsz(4)*2/3;
        gap = 5;
        h.fig = uifigure('Position',round([(scsz(3)-fwidth)/2, scsz(4)/6, fwidth, fheight]),...
            'Name','Select DICOMcatalog data for processing:');%,'CloseRequestFcn',@cancel_callback);
        h.tabgroup = uitabgroup(h.fig,'Position',[0,0,fwidth,fheight]);
        h.tab2 = uitab(h.tabgroup,'Title','Options');
        h.tab1 = uitab(h.tabgroup,'Title','Scans');

        %% Tab 1: tables for selecting images, and execution buttons
        h.table_select = uitable(h.tab1,'Position',[ 1 , 20 , fwidth , fheight-130 ],'CellEditCallback',@editCell);
        h.table_filter = uitable(h.tab1,'Position',[ 1 , fheight-115 , fwidth , 90 ],'CellEditCallback',@applyFilter);

        uibutton(  h.tab1,'Position',[1,1,50,20],'Text','Clear','BackgroundColor','blue','ButtonPushedFcn',@clearTags);
        
        uilabel(h.tab1,'Position',[fwidth/2-320,1,135,20],'HorizontalAlignment','right','Text','Groups per page: ');
        h.edit_page_groups = uieditfield(h.tab1,'numeric','Position',[fwidth/2-185,1,50,20],...
            'Editable',1,'Value',gp_per_page,'HorizontalAlignment','center','ValueChangedFcn',@setPage,'Tag','edit_gp_per');
        h.text_groupN = uilabel(h.tab1,'Position',[fwidth/2-135,1,60,20],'Text',' of 0');
        
        uilabel(h.tab1,'Position',[fwidth/2-105,1,60,20],'HorizontalAlignment','right','Text','Page: ');
        uibutton(  h.tab1,'Position',[fwidth/2-45,1,20,20],'Text','<','Tag','decr','ButtonPushedFcn',@setPage);
        h.text_page = uieditfield(h.tab1,'numeric','Position',[fwidth/2-25,1,50,20],...
            'Editable',1,'Value',0,'HorizontalAlignment','center','ValueChangedFcn',@setPage,'Tag','text_page');
        uibutton(  h.tab1,'Position',[fwidth/2+25,1,20,20],'Text','>','Tag','incr','ButtonPushedFcn',@setPage);
        h.text_pageN = uilabel(h.tab1,'Position',[fwidth/2+45,1,60,20],'Text',' of 0');
        
        uibutton(  h.tab1,'Position',[fwidth-50,1,50,20],...
            'Text','Done','BackgroundColor','green','ButtonPushedFcn',@done_callback);
        uibutton(  h.tab1,'Position',[fwidth-100-gap,1,50,20],...
            'Text','Cancel','BackgroundColor','red','ButtonPushedFcn',@cancel_callback);

        %% Tab 2: Settings and options
        uitextarea(h.tab2,'Position',[5,   fheight-50, 100,        20],'Editable',0,'Value','uniquename:');
        h.edit_username = uieditfield(h.tab2,'Position',[110, fheight-50, 300, 20],'Editable',1,...
            'ValueChangedFcn',@setOpts,'Tag','username','Value',opts.username);
        uibutton(h.tab2,  'Position',[5,   fheight-75, 100,        20],'Text','Select Catalog:','ButtonPushedFcn',@selectCatalog);
        h.text_cat = uitextarea(h.tab2,'Position',[110, fheight-75, fwidth-115, 20],'Editable',0);
        uibutton(h.tab2,  'Position',[5,   fheight-100, 100,        20],'Text','Save To:','ButtonPushedFcn',@setSavePath);
        h.text_save = uitextarea(h.tab2,'Position',[110, fheight-100, fwidth-115, 20],'Editable',0);

        h.panel_modules = uipanel(h.tab2,'Position',[5, fheight-330, 215, 225],'Title','Modules');
        uibutton(h.panel_modules,  'Position',[5,   185, 100, 20],'Text','Select All','ButtonPushedFcn',@selectAll);
        uibutton(h.panel_modules,  'Position',[110, 185, 100, 20],'Text','Clear All','ButtonPushedFcn',@clearAll);
        h.unreg =   uicheckbox(h.panel_modules,'Position',[5,   160, 200, 20],'Text','Unreg',...
            'Value',opts.unreg,'ValueChangedFcn',@setOpts,'Tag','unreg');
        h.airway =  uicheckbox(h.panel_modules,'Position',[5,   135, 200, 20],'Text','Airways',...
            'Value',opts.airway,'ValueChangedFcn',@setOpts,'Tag','airway');
        h.scatnet = uicheckbox(h.panel_modules,'Position',[5,   110, 200, 20],'Text','ScatNet',...
            'Value',opts.scatnet,'ValueChangedFcn',@setOpts,'Tag','scatnet');
        h.vessel =  uicheckbox(h.panel_modules,'Position',[5,    85, 200, 20],'Text','Vessels',...
            'Value',opts.vessel,'ValueChangedFcn',@setOpts,'Tag','vessel');
        h.reg =     uicheckbox(h.panel_modules,'Position',[5,    60, 200, 20],'Text','Registration',...
            'Value',opts.reg,'ValueChangedFcn',@setOpts,'Tag','reg');
        h.prm =     uicheckbox(h.panel_modules,'Position',[5,    35, 200, 20],'Text','PRM',...
            'Value',opts.prm,'ValueChangedFcn',@setOpts,'Tag','prm');
        h.tprm =    uicheckbox(h.panel_modules,'Position',[5,    10, 200, 20],'Text','tPRM',...
            'Value',opts.tprm,'ValueChangedFcn',@setOpts,'Tag','tprm');

        h.panel_reg = uipanel(h.tab2,'Position',[225, fheight-330, 215, 225],'Title','Reg Options');
        h.reg_seg =  uicheckbox(h.panel_reg,'Position',[5, 185, 200, 20],'Text','Tranform Segmentation',...
            'Value',opts.reg_seg,'ValueChangedFcn',@setOpts,'Tag','reg_seg','Enable','off');
        h.quickreg = uicheckbox(h.panel_reg,'Position',[5,160,200,20],'Text','Quick Registration',...
            'Value',opts.quickreg,'ValueChangedFcn',@setOpts,'Tag','quickreg');
        h.jac =      uicheckbox(h.panel_reg,'Position',[5,135,200,20],'Text','|Jacobian|',...
            'Value',opts.jac,'ValueChangedFcn',@setOpts,'Tag','jac');
        h.jacmat =   uicheckbox(h.panel_reg,'Position',[5,110,200,20],'Text','Jacobian Matrix',...
            'Value',opts.jacmat,'ValueChangedFcn',@setOpts,'Tag','jacmat');
        h.def =      uicheckbox(h.panel_reg,'Position',[5,85,200,20],'Text','Deformation Fields',...
            'Value',opts.def,'ValueChangedFcn',@setOpts,'Tag','def');
        h.Tvoi =  uicheckbox(h.panel_reg,'Position',[5,60,200,20],'Text','Transform VOI',...
            'Value',opts.Tvoi,'ValueChangedFcn',@setOpts,'Tag','Tvoi');
        h.Tsurf = uicheckbox(h.panel_reg,'Position',[5,35,200,20],'Text','Transform Surface',...
            'Value',opts.Tsurf,'ValueChangedFcn',@setOpts,'Tag','Tsurf');
        
        h.panel_yacta = uipanel(h.tab2,'Position',[445, fheight-330, 215, 225],'Title','YACTA Options');
%         h.y_airways
%         h.y_parenchyma
        
        h.panel_other = uipanel(h.tab2,'Position',[665, fheight-330, 215, 225],'Title','Other Options');
        h.orient_check = uicheckbox(h.panel_other,'Position',[5,185, 200, 20],'Text','Auto-correct Orientation',...
            'Value',opts.orient_check,'ValueChangedFcn',@setOpts,'Tag','orient_check');
        uilabel(h.panel_other,'Position',[5, 160, 70, 20],'HorizontalAlignment','right','Text','Parallel Pool Size:');
        h.edit_npar = uieditfield(h.panel_other,'numeric','Position',[75,160,50,20],...
            'Editable',1,'Value',opts.par_size,'HorizontalAlignment','center','ValueChangedFcn',@setOpts,'Tag','par_size');
        
    end
    function setPage(hObject,eventData)
        if isnumeric(hObject)
            hObject = round(hObject);
            if hObject>0 && hObject<pageN
                curr_page = hObject;
            else
                warning('Invalid numerical input. Max page = %s',pageN);
            end
        elseif isprop(hObject,'Tag') % incr/decr button
            switch hObject.Tag
                case 'incr'
                    curr_page = curr_page + 1;
                case 'decr'
                    curr_page = curr_page - 1;
                case 'text_page'
                    val = round(eventData.Value);
                    if val>0 && val<=ngroups
                        curr_page = val;
                    else
                        h.text_page.Value = eventData.PreviousValue;
                        return;
                    end
                case 'edit_gp_per'
                    val = round(eventData.Value);
                    if val>0
                        gp_per_page = val;
                        pageN = ceil(ngroups/gp_per_page);
                        h.edit_gp_per_page.Value = gp_per_page;
                        h.text_pageN.Text = sprintf('of %u',pageN);
                    else
                        h.edit_gp_per_page.Value = eventData.PreviousValue;
                        return;
                    end
                otherwise
                    return;
            end
        else % Invalid input
            return;
        end
        h.text_page.Value = curr_page;
        page_ic = ismember(ugroups_ic,(gp_per_page*(curr_page-1))+1:min(ngroups,gp_per_page*curr_page));
        h.table_select.Data = C(page_ic,:);
        checkValid();
    end
    function applyFilter(~,eventdata)
        nscans = size(C,1);
        
        if nargin==0 % filter all
            rowi = 1:2;
        elseif nargin==2 % filter single
            rowi = eventdata.Indices(1);
        end
        
        % Set tags:
        for i = rowi
            ei_str = h.table_filter.Data{i,1}{1};
            if any(~cellfun(@isempty,table2cell(h.table_filter.Data(i,2:end))))
                TF = true(nscans,1);
                filtnames = h.table_filter.ColumnName(2:end);
                filtnames(cellfun(@isempty,table2cell(h.table_filter.Data(i,2:end)))) = [];
                for j = 1:numel(filtnames)
                    filtval = h.table_filter.Data.(filtnames{j}){i};
                    if iscell(C.(filtnames{j}))
                        % String comparison
                        TF = TF & contains(C.(filtnames{j}),filtval);
                    elseif isnumeric(C.(filtnames{j}))
                        % Numeric comparison
                        TF = TF & (C.(filtnames{j})==str2double(filtval));
                    end
                end
            else
                TF = false(nscans,1);
            end
            C.(ei_str) = TF;
            h.table_select.Data.(ei_str) = TF(page_ic);
        end

        checkValid;
    end

    function editCell(hObject,eventdata)
        if ismember(eventdata.Indices(2),[1,2])
            % Edited EXP/INS Tag
            C{eventdata.Indices(1)+find(page_ic,1)-1,eventdata.Indices(2)} = eventdata.NewData;
            checkValid;
        elseif isempty(regexp(eventdata.NewData , '[/\*:?"<>|]', 'once')) % Edited UMlabel
            % Set ID for all scans in this case group:
            page_group_ic = ugroups_ic(page_ic);
            sub_ic = page_group_ic==page_group_ic(eventdata.Indices(1));
            hObject.Data(sub_ic,3) = {eventdata.NewData};
            C(page_ic(sub_ic),3) = {eventdata.NewData};
        else
            warning('Invalid PatientID, try again.');
            hObject.Data{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.PreviousData;
        end
    end

    function checkValid(~)
        page_group_ic = ugroups_ic(page_ic);
        ugroup = unique(page_group_ic);
        for i = 1:numel(ugroup)
            irow = find(page_group_ic==ugroup(i));
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
        C(page_ic,1:2) = {false};
        checkValid;
    end
    function done_callback(~,~)
        % Check that each timepoint has Ins/Exp selected (and only one of each)
        if any(gp_valid>2)
            warning('No case can have more than one Exp/Ins selected.');
        else
            % Set up options structure
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
        str = hObject.Tag;
        if ismember(str,fieldnames(opts))
            opts.(str) = hObject.Value;
        end
    end
    function str = selectCatalog(str,~)
        newC = [];
        if nargin==0 || isempty(str) || ~ischar(str)
            str = uigetdir('Select folder for processing.');
            if isnumeric(str) % Cancel button was selected
                return;
            end
        end
        if istable(str)
            newC = str;
            str = '';
        elseif ischar(str) && isfolder(str)
            % Find existing DICOMcatalog files:
            fnames = dir(fullfile(str,'*DICOMcatalog*.csv'));
            answer = listdlg('PromptString','Select a catalog:',...
                'ListString',{'* Create New *',fnames.name});
            if isempty(answer)
                str = '';
            elseif answer == 1
                % Generate DICOM catalog:
                newC = dcmCatalog(str);
                str = fullfile(str,'DICOMcatalog.csv');
            elseif answer > 1
                str = fullfile(fnames(answer-1).folder,fnames(answer-1).name);
            end
        end
        if isempty(newC) && ischar(str) && exist(str,'file')==2
            iopt = detectImportOptions(str);
            newC = readtable(str,iopt);
        end
        if ~isempty(newC) % Set up the tables and data:
            %% Validate table input:
            colnames = newC.Properties.VariableNames;
            i_member = ismember(req_fields,colnames);
            if all(ismember(req_fields,colnames))
                C = newC;
            else
                error(['Missing required fields: ',strjoin(req_fields(~i_member),', ')]);
            end
            
            %% Set up display data:
            % Re-order columns for easier reading:
            C = movevars(C,req_fields,'Before',1);

            % Add Tag columns:
            nscans = size(C,1);
            C = addvars(C,false(nscans,1),false(nscans,1),'Before',1,'NewVariableNames',{'Exp','Ins'});
            colnames = C.Properties.VariableNames;
            
            % Check for saved selections:
            if ismember('Tag',colnames)
                C.Exp(strcmp(C.Tag,'Exp')) = true;
                C.Ins(strcmp(C.Tag,'Ins')) = true;
                C = removevars(C,'Tag');
            end
            if ismember('UMlabel',colnames)
                C = movevars(C,'UMlabel','After','Ins');
            else
                UMlabel = C.PatientName;
                ind = cellfun(@(x)isempty(x)||strcmp(x,'Anonymous'),UMlabel);
                UMlabel(ind) = C.PatientID(ind);
                C = addvars(C,UMlabel,'After','Ins');
            end

            % Remove empty data and sort array
            % Find groups of scans with unique PatientName and StudyDate
            if isnumeric(C.StudyDate)
                C.StudyDate = cellfun(@num2str,num2cell(C.StudyDate),'UniformOutput',false);
            end
            if ismember('CaseNumber',colnames)
                C = removevars(C,'CaseNumber');
            end
            % Remove cases with no identifiers
            C(cellfun(@isempty,C.UMlabel),:) = [];
            % Sort by identifiers
            C = sortrows(C,{'StudyDate','StudyID','PatientName','SeriesNumber'});
            % Find case groupings
            if isnumeric(C.PatientName)
                C.PatientName = cellfun(@num2str,num2cell(C.PatientName),'UniformOutput',false);
            end
            [ugroups,~,ugroups_ic] = unique(strcat(C.StudyDate,C.StudyID,C.PatientName,C.StudyID));
            ngroups = numel(ugroups);

            gp_valid = zeros(ngroups,1); % Initialize group validation
        
            % Set page info
            pageN = ceil(ngroups/gp_per_page);
            h.text_pageN.Text = sprintf('of %d',pageN);
            h.text_groupN.Text = sprintf('of %d',ngroups);
            
            %% Set filter table:
            nfields = numel(C.Properties.VariableNames);
            colEdit = false(1,nfields); colEdit([1,2,3]) = true;
            colWidth = repmat({'auto'},1,nfields); colWidth(1:2) = {40,40};
            nv = size(C,2)-3;
            fC = table('Size',[2,nv],...
                'VariableTypes',repmat({'cellstr'},1,nv),...
                'VariableNames',C.Properties.VariableNames(4:end));
            fC = addvars(fC,{'Exp';'Ins'},'Before',1,'NewVariableNames',{'Tag'});
            set(h.table_filter,'Data',fC,'ColumnEditable',[false,true(1,nv)],'ColumnWidth',[{80},colWidth(3:end)]);
            
            %% Set table data:
            setPage(1);
            
            % Set table properties
            set(h.table_select,'ColumnEditable',colEdit,'ColumnWidth',colWidth);
            
            pause(1);
            checkValid;
            
        end
        opts.dcm_path = str;
        h.text_cat.Value = str;
        figure(h.fig);
    end
    function setSavePath(tpath,~)
        if nargin==0 || isempty(tpath) || ~ischar(tpath)
            tpath = uigetdir('Select folder for saving results:');
        end
        if ~isempty(tpath)
            opts.save_path = tpath;
            h.text_save.Value = tpath;
        end
        figure(h.fig);
    end
    function selectAll(~,~)
        set([h.unreg,h.airway,h.scatnet,h.vessel,h.reg,h.prm,h.tprm],'Value',1);
    end
    function clearAll(~,~)
        set([h.unreg,h.airway,h.scatnet,h.vessel,h.reg,h.prm,h.tprm],'Value',0);
    end
end
