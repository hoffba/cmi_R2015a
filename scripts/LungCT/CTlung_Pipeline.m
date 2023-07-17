function CTlung_Pipeline(varargin)
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
opts = init_pipeline_opts;
          
C = [];          % full catalog table

% Scan groups
ugroups_ic = []; % unique group indices
ngroups = 0;     % number of cases available with unique ID and date
gp_valid = [];   
gp_caselabel = {};  % Full label for saving case files - <UMlabel>_<StudyDate>
gp_sel = [];        % TF - which groups are selected

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

% Preview selection
preview_ind = [];
cmiObj = [];

% Help figure
h_help = [];


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

tagstr = {' ','Exp','Ins','ExpLabel','InsLabel'};

%% Initialize the GUI:
h = initFig;
drawnow;

if ~isempty(opts.dcm_path)
    selectCatalog(opts.dcm_path);
end
if ~isempty(opts.save_path)
    setSavePath(opts.save_path);
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
        h.fig = uifigure('Position',round([(scsz(3)-fwidth)/2, scsz(4)/6, fwidth, fheight]),...
            'Name','Select DICOMcatalog data for processing:');
        h.tabgroup = uitabgroup(h.fig,'Position',[0,0,fwidth,fheight]);
        h.tab_opts = uitab(h.tabgroup,'Title','Options');
        h.tab_scans = uitab(h.tabgroup,'Title','Scans');

    %% Tab 1: Settings and options
        uilabel(h.tab_opts,   'Position',[5,   fheight-50, 100, 20],'Text','Uniquename:','FontWeight','bold');
        h.edit_username = uieditfield(h.tab_opts,'Position',[110, fheight-50, 300, 20],'Editable',1,...
            'ValueChangedFcn',@setOpts,'Tag','username','Value',opts.username);
        uibutton(h.tab_opts,  'Position',[fwidth-45,  fheight-50, 40, 20],'Text','Help','ButtonPushedFcn',@pipelineHELP,'BackgroundColor',[.73,.92,1]);
        uilabel(h.tab_opts,   'Position',[5,   fheight-75, 70, 20],'Text','Catalog: ','FontWeight','bold');
        uibutton(h.tab_opts,  'Position',[80,  fheight-75, 50, 20],'Text','Select','ButtonPushedFcn',@selectCatalog,'Tag','cat_select');
        uibutton(h.tab_opts,  'Position',[135, fheight-75, 50, 20],'Text','New','ButtonPushedFcn',@selectCatalog,'Tag','cat_new');
        uibutton(h.tab_opts,  'Position',[190, fheight-75, 50, 20],'Text','Save','ButtonPushedFcn',@saveCatalog);
        h.text_cat = uitextarea(h.tab_opts,'Position',[80, fheight-100, fwidth-85, 20],'Editable',0);

        uilabel(h.tab_opts,   'Position',[5,   fheight-125,  70, 20],'Text','Results:','FontWeight','bold');
        uibutton(h.tab_opts,  'Position',[80,   fheight-125, 50,20],'Text','Select','ButtonPushedFcn',@setSavePath);
        h.text_save = uitextarea(h.tab_opts,'Position',[80, fheight-150, fwidth-85, 20],'Editable',0);

    % Panel: Modules
        h.panel_modules = uipanel(h.tab_opts,'Position',[5, fheight-440, 215, 285],'Title','Modules','FontWeight','bold');
        uibutton(h.panel_modules,  'Position',[5,   245, 100, 20],'Text','Select All','ButtonPushedFcn',@selectAll);
        uibutton(h.panel_modules,  'Position',[110, 245, 100, 20],'Text','Clear All','ButtonPushedFcn',@clearAll);
        h.unreg =   uicheckbox(h.panel_modules,'Position',[5,   210, 200, 20],'Text','Unreg',...
            'Value',opts.unreg,'ValueChangedFcn',@setOpts,'Tag','unreg');
        h.airway =  uicheckbox(h.panel_modules,'Position',[5,   185, 200, 20],'Text','Airways',...
            'Value',opts.airway,'ValueChangedFcn',@setOpts,'Tag','airway');
        h.scatnetAT = uicheckbox(h.panel_modules,'Position',[5,   160, 200, 20],'Text','ScatNet-AT',...
            'Value',opts.scatnetAT,'ValueChangedFcn',@setOpts,'Tag','scatnetAT');
        h.scatnetAT_PEDS = uicheckbox(h.panel_modules,'Position',[5,   135, 200, 20],'Text','ScatNet-AT-Peds',...
            'Value',opts.scatnetAT_PEDS,'ValueChangedFcn',@setOpts,'Tag','scatnetAT_PEDS');
        h.vessel =  uicheckbox(h.panel_modules,'Position',[5,    110, 200, 20],'Text','Vessels',...
            'Value',opts.vessel,'ValueChangedFcn',@setOpts,'Tag','vessel');
        h.reg =     uicheckbox(h.panel_modules,'Position',[5,    85, 200, 20],'Text','Registration',...
            'Value',opts.reg,'ValueChangedFcn',@setOpts,'Tag','reg');
        h.scatnetEmph = uicheckbox(h.panel_modules,'Position',[5,   60, 200, 20],'Text','ScatNet-Emph',...
            'Value',opts.scatnetEmph,'ValueChangedFcn',@setOpts,'Tag','scatnetEmph');
        h.prm =     uicheckbox(h.panel_modules,'Position',[5,    35, 200, 20],'Text','PRM',...
            'Value',opts.prm,'ValueChangedFcn',@setOpts,'Tag','prm');
        h.tprm =    uicheckbox(h.panel_modules,'Position',[5,    10, 200, 20],'Text','tPRM',...
            'Value',opts.tprm,'ValueChangedFcn',@setOpts,'Tag','tprm');

    % Panel: Registration Options
        h.panel_reg = uipanel(h.tab_opts,'Position',[225, fheight-415, 215, 260],'Title','Reg Options','FontWeight','bold');
        h.dBlood =  uicheckbox(h.panel_reg,'Position',[5, 185, 200, 20],'Text','Blood density change',...
            'Value',opts.dBlood,'ValueChangedFcn',@setOpts,'Tag','dBlood');
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

    % Panel: Other Options
        h.panel_other = uipanel(h.tab_opts,'Position',[445, fheight-415, 215, 260],'Title','Other Options','FontWeight','bold');
        h.orient_check = uicheckbox(h.panel_other,'Position',[5,185, 200, 20],'Text','Auto-correct Orientation',...
            'Value',opts.orient_check,'ValueChangedFcn',@setOpts,'Tag','orient_check');
        h.peds_check = uicheckbox(h.panel_other,'Position',[5,   160, 200, 20],'Text','Pediatrics',...
            'Value',opts.peds,'ValueChangedFcn',@setOpts,'Tag','peds');
        
    % Panel: Cluster Options
        h.panel_cluster = uipanel(h.tab_opts,'Position',[665, fheight-415, 215, 260],'Title','Cluster Options','FontWeight','bold');
        bg = uibuttongroup(h.panel_cluster,'Position',[5,110,205,95],'Title','Cluster:',...
            'SelectionChangedFcn',@setOpts,'Tag','cluster');
        h.radio_GL = uiradiobutton(bg,'Position',[5,55,195,20],'Text','Great Lakes','Tag','GL');
        h.radio_batch = uiradiobutton(bg,'Position',[5,30,195,20],'Text','Local Batch','Tag','batch');
        h.radio_debug = uiradiobutton(bg,'Position',[5,5,195,20],'Text','Debug','Tag','debug');
        h.(sprintf('radio_%s',opts.cluster)).Value = true;
        uilabel(h.panel_cluster,'Position',[5, 80, 120, 20],'Text','Parallel Pool Size:');
        h.edit_npar = uieditfield(h.panel_cluster,'numeric','Position',[125,80,50,20],'Enable',false,...
            'Editable',1,'Value',opts.par_size,'HorizontalAlignment','center','ValueChangedFcn',@setOpts,'Tag','par_size');
        uilabel(h.panel_cluster,'Position',[5, 55, 120, 20],'Text','Partition Type:');
        h.partition = uidropdown(h.panel_cluster,'Position',[125,55,80,20],...
            'Items',{'auto','standard','largemem'},'ValueChangedFcn',@setOpts,'Tag','partition');
        uilabel(h.panel_cluster,'Position',[5, 30, 120, 20],'Text','Memory / process: ');
        h.edit_mem = uieditfield(h.panel_cluster,'numeric','Position',[125, 30, 50, 20],...
            'Editable',1,'Value',opts.mem,'HorizontalAlignment','center','ValueChangedFcn',@setOpts,'Tag','mem');
        uilabel(h.panel_cluster,'Position',[5, 5, 120, 20],'Text','Number of Nodes: ');
        h.nnodes = uidropdown(h.panel_cluster,'Position',[125, 5, 80, 20],...
            'Items',{'1','2','3','4','5'},'ItemsData',1:5,'ValueChangedFcn',@setOpts,'Tag','nnodes');
                
    %% Tab 2: tables for selecting images, and execution buttons
        % Filter selection:
        uilabel(h.tab_scans,'Position',[5,fheight-50,50,20],'Text','Filter:');
        h.filt_tag_dd = uidropdown(h.tab_scans,'Position',[60,fheight-50,100,20],'Items',tagstr); % Select Tag to filter for
        h.filt_fld_dd = uidropdown(h.tab_scans,'Position',[165,fheight-50,200,20],'Items',{});
        uibutton(h.tab_scans,'Position',[370,fheight-50,50,20],'Text','Apply','ButtonPushedFcn',@applyFilter);
        h.filt_bg = uibuttongroup(h.tab_scans,'Position',[425,fheight-50,135,20],'Tag','filtselect');
        h.radio_filt_all = uiradiobutton(h.filt_bg,'Position',[5,0,60,20],'Text','All','Tag','filt_all');
        h.radio_filt_page = uiradiobutton(h.filt_bg,'Position',[70,0,60,20],'Text','Page','Tag','filt_page');
        h.edit_filter = uieditfield(h.tab_scans,'Position',[60,fheight-75,500,20]);

        % Scans table
        h.table_select = uitable(h.tab_scans,'Position',[ 1 , 20 , fwidth , fheight-100 ],'CellEditCallback',@editCell,'CellSelectionCallback',@selectCell);
        uibutton(  h.tab_scans,'Position',[1,1,50,20],...
            'FontWeight','bold','Text','Clear','BackgroundColor',[.73,.92,1],'ButtonPushedFcn',@clearTags);
        
        uilabel(   h.tab_scans,'Position',[fwidth/2-320,1,135,20],'HorizontalAlignment','right','Text','Groups per page: ');
        h.edit_page_groups = uieditfield(h.tab_scans,'numeric','Position',[fwidth/2-185,1,50,20],...
            'Editable',1,'Value',gp_per_page,'HorizontalAlignment','center','ValueChangedFcn',@setPage,'Tag','edit_gp_per');
        h.text_groupN = uilabel(h.tab_scans,'Position',[fwidth/2-135,1,60,20],'Text',' of 0');
        uilabel(   h.tab_scans,'Position',[fwidth/2-105,1,60,20],'HorizontalAlignment','right','Text','Page: ');
        uibutton(  h.tab_scans,'Position',[fwidth/2-45,1,20,20],'Text','<','Tag','decr','ButtonPushedFcn',@setPage);
        h.text_page = uieditfield(h.tab_scans,'numeric','Position',[fwidth/2-25,1,50,20],...
            'Editable',1,'Value',0,'HorizontalAlignment','center','ValueChangedFcn',@setPage,'Tag','text_page');
        uibutton(  h.tab_scans,'Position',[fwidth/2+25,1,20,20],'Text','>','Tag','incr','ButtonPushedFcn',@setPage);
        h.text_pageN = uilabel(h.tab_scans,'Position',[fwidth/2+45,1,60,20],'Text',' of 0');
        
        uibutton(  h.tab_scans,'Position',[fwidth/2+130,1,75,20],...
            'Text','Preview','BackgroundColor','yellow','ButtonPushedFcn',@preview);
        
        uibutton(  h.tab_scans,'Position',[fwidth-50,1,50,20],...
            'Text','Start','BackgroundColor','green','FontWeight','bold','ButtonPushedFcn',@run);


    end
    function setPage(hObject,eventData)
        if isnumeric(hObject)
            hObject = round(hObject);
            if hObject>0 && hObject<=pageN
                newpage = hObject;
            else
                warning('Invalid numerical input. Max page = %s',pageN);
            end
        elseif isprop(hObject,'Tag') % incr/decr button
            switch hObject.Tag
                case 'incr'
                    newpage = min(curr_page + 1,pageN);
                case 'decr'
                    newpage = max(curr_page - 1,1);
                case 'text_page'
                    val = round(eventData.Value);
                    if val>0 && val<=pageN
                        newpage = val;
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
        if newpage ~= curr_page
            curr_page = newpage;
            h.text_page.Value = curr_page;
            page_ic = ismember(ugroups_ic,(gp_per_page*(curr_page-1))+1:min(ngroups,gp_per_page*curr_page));
            h.table_select.Data = table2cell(C(page_ic,:));
            h.table_select.ColumnName = C.Properties.VariableNames;
            nv = size(C,2);
            set(h.table_select,'ColumnEditable',[true(1,2),false(1,nv-2)],...
                               'ColumnFormat',{tagstr}); % Tag options
            checkValid();
        end
    end
    function applyFilter(~,~)
        filt_tag = h.filt_tag_dd.Value;
        filt_fld = h.filt_fld_dd.Value;
        page_flag = strcmp(h.filt_bg.SelectedObject.Text,'Page');
        filt_value = h.edit_filter.Value;

        % Extract existing Tag and selected field
        tC = C.(filt_fld);
        tags = C.Tag;
        if page_flag
            tC = tC(page_ic);
            tags = tags(page_ic);
        end

        % Apply filter
        TF = false(numel(tC),1);
        if iscell(tC)
            TF = ~cellfun(@isempty,regexp(tC,filt_value));
        elseif isnumeric(tC)
            TF = tC==str2double(filt_value);
        end

        if any(TF)
            TF_tag = strcmp(tags,filt_tag);
            if any(TF_tag)
                % Apply filter to exiting selections
                tags(TF_tag) = {''};
                TF = TF & TF_tag;
            end
            tags(TF) = {filt_tag};
            
            % Set Tags in table
            if page_flag
                C.Tag(page_ic) = tags;
            else
                C.Tag = tags;
            end

            % Update current page
            h.table_select.Data(:,1) = tags;

            % Validate new selections
            checkValid;
        end
    end

    function editCell(hObject,eventdata)
        C_ind = eventdata.Indices(1)+find(page_ic,1)-1;
        gp = ugroups_ic(C_ind);
        if eventdata.Indices(2)==1
            % Edited Tag
            newval = eventdata.NewData;
            if strcmp(newval,' ')
                newval = '';
            end
            C.Tag{C_ind,eventdata.Indices(2)} = newval;
            
            % Update case selection status
            gp_sel(gp) = any(ismember(C.Tag(ugroups_ic==gp),{'Exp','Ins'}));
            
        elseif isempty(regexp(eventdata.NewData , '[/\*:?"<>|]', 'once')) 
            % Edited UMlabel
            % Set ID for all scans in this case group:
            str = regexprep(eventdata.NewData,' ','_'); % replace spaces with underscore
            hObject.Data(ugroups_ic(page_ic)==gp,2) = {str}; % update page table for case
            C.UMlabel(ugroups_ic==gp) = {eventdata.NewData}; % update C data for case
            
            % Update case name
            gp_caselabel{gp} = strcat(str,'_',C.StudyDate{C_ind});
            
        else
            warning('Invalid PatientID, try again.');
            hObject.Data{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.PreviousData;
        end
        checkValid;
    end

    function checkValid(~)
        page_group_ic = ugroups_ic(page_ic);
        ugroup = unique(page_group_ic);
        for i = 1:numel(ugroup)
            gp_ind = i+(curr_page-1)*gp_per_page;
            irow = find(page_group_ic==ugroup(i));
            n_exp = nnz(strcmp(h.table_select.Data(irow,1),'Exp'));
            n_ins = nnz(strcmp(h.table_select.Data(irow,1),'Ins'));
            gp_valid(gp_ind) = n_exp + n_ins;
            if gp_valid(gp_ind)
                % Check for multiple selections:
                if (n_exp>1) || (n_ins>1)
                    gp_valid(gp_ind) = 3;
                end
                % Check for duplicate UMlabel: gp_umlabel
                if nnz(strcmp(gp_caselabel(gp_ind),gp_caselabel(gp_sel)))>1
                    gp_valid(gp_ind) = 4;
                end
            end
            switch gp_valid(gp_ind)
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
        h.table_select.Data(:,1) = {''};
        C.Tag(page_ic) = {''};
        checkValid;
    end
    function setOpts(hObject,eData)
        tag = hObject.Tag;
        switch tag
            case 'cluster'
                h.edit_npar.Enable = strcmp(eData.NewValue.Tag,'batch');
                set([h.partition,h.edit_mem,h.nnodes],'Enable',strcmp(eData.NewValue.Tag,'GL'));
                opts.cluster = eData.NewValue.Tag;
            case 'partition'
                Ncases = numel(unique(ugroups_ic(C.Exp|C.Ins)));
                stat = checkSLURM(numcases,eData.Value,opts.mem,opts.nnodes);
                opts.partition = eData.Value;
            case 'mem'
                opts.mem = eData.Value;
            case 'nnodes'
                opts.nnodes = eData.Value;
            otherwise
                if ismember(tag,fieldnames(opts))
                    opts.(tag) = hObject.Value;
                end
        end
    end
    function fn = selectCatalog(fn,~)
        if isa(fn,'matlab.ui.control.Button')
            switch fn.Tag
                case 'cat_select'
                    [fn,path] = uigetfile('*.csv','Select Catalog File:');
                    if ischar(fn)
                        fn = fullfile(path,fn);
                    else
                        return;
                    end
                case 'cat_new'
                    fn = uigetdir(fileparts(opts.dcm_path),'Select folder for processing.');
                otherwise
                    return;
            end
        elseif isempty(fn) || ~ischar(fn)
            return;
        end
        if istable(fn)
            newC = fn;
            fn = '';
        elseif ischar(fn) && isfolder(fn)
            % If a folder was input, generate new catalog
            newC = catalog_data(fn);
            fn = fullfile(fn,'Pipeline_catalog.csv');
        elseif ischar(fn) && exist(fn,'file')==2
            iopt = detectImportOptions(fn);
            newC = readtable(fn,iopt);
        else
            return;
        end
        if ~isempty(newC) % Set up the tables and data:
            
            %% Validate table input:
            if ismember('Directory',newC.Properties.VariableNames)
                newC = renamevars(newC,'Directory','DataPath');
            end
            i_member = ismember(req_fields,newC.Properties.VariableNames);
            if all(i_member)
                C = newC;
            else
                error(['Missing required fields: ',strjoin(req_fields(~i_member),', ')]);
            end
            
            %% Set up display data:
            % Make sure Tag column is first
            if ismember('Tag',C.Properties.VariableNames)
                C = movevars(C,'Tag','Before',1);
            else
                C = addvars(C,cellstr(size(C,1),1),'Before',1,'NewVariableNames',{'Tag'});
            end

            % Re-order columns for easier reading:
            C = movevars(C,req_fields,'After',1);
            
            if ismember('UMlabel',C.Properties.VariableNames)
                C = movevars(C,'UMlabel','After',1);
            else
                UMlabel = C.PatientName;
                
                % Try setting empty/anonymouse names using PatientID instead
                ind = cellfun(@(x)isempty(x)||strcmp(x,'Anonymous'),UMlabel);
                UMlabel(ind) = C.PatientID(ind);
                
                C = addvars(C,UMlabel,'After',1);
            end
            if isnumeric(C.UMlabel)
                C.UMlabel = arrayfun(@num2str,C.UMlabel,'UniformOutput',false);
            end

            % Remove empty data and sort array
            % Find groups of scans with unique PatientName and StudyDate

            % Force these fields to be cellstr
            tlabel = {'UMlabel','StudyID','StudyDate','PatientName'};
            for i = 1:numel(tlabel)
                if ismember(tlabel{i},C.Properties.VariableNames) && isnumeric(C.(tlabel{i}))
                    C.(tlabel{i}) = arrayfun(@num2str,C.(tlabel{i}),'UniformOutput',false);
                end
            end
            
            % Remove spaces from UMlabel
            C.UMlabel = cellfun(@(x)regexprep(x,' ','_'),C.UMlabel,'UniformOutput',false);
            
            if ismember('CaseNumber',C.Properties.VariableNames)
                C = sortrows(C,{'StudyDate','StudyID','PatientName','CaseNumber','SeriesNumber'});
                [~,~,ugroups_ic] = unique(C.CaseNumber);        % Find unique CaseNumbers
                ugroups_ic = [0;cumsum(diff(ugroups_ic)~=0)]+1; % Re-number case groups
                C.CaseNumber = ugroups_ic;
                [~,ugroups_ia] = unique(ugroups_ic);            % Find new case group locations
            else
                % Remove cases with no identifiers
                C(cellfun(@isempty,C.UMlabel),:) = [];
                % Sort by identifiers
                C = sortrows(C,{'StudyDate','StudyID','PatientName','SeriesNumber'});
                % Find case groupings
                [~,ugroups_ia,ugroups_ic] = unique(strcat(C.StudyDate,C.StudyID,C.PatientName));
            end
            ngroups = numel(ugroups_ia);

            % Initialize group info
            gp_valid = zeros(ngroups,1); % group validation
            gp_caselabel = strcat(C.UMlabel(ugroups_ia),'_',C.StudyDate(ugroups_ia));
            gp_sel = false(ngroups,1);
            for i = 1:ngroups
                gp_sel(i) = any(ismember(C.Tag(ugroups_ic==i),{'Exp','Ins'}));
            end
            
            % Fix DICOM location
            if ischar(fn)
                cat_path = fileparts(fn);
                for ifn = 1:size(C,1)
                    [~,cat_folder] = fileparts(cat_path);
                    tpath = extractAfter(C.DataPath{ifn},cat_folder);
                    C.DataPath{ifn} = fullfile(cat_path,tpath);
                end
            end
        
            % Set page info
            pageN = ceil(ngroups/gp_per_page);
            h.text_pageN.Text = sprintf('of %d',pageN);
            h.text_groupN.Text = sprintf('of %d',ngroups);
            
            %% Set filter table:
            set(h.filt_fld_dd,'Items',C.Properties.VariableNames(3:end))
            
            %% Set table data:
            setPage(1);
            
            pause(1);
            checkValid;
            
        end
        opts.dcm_path = fn;
        h.text_cat.Value = fn;
        figure(h.fig);
    end
    function setSavePath(tpath,~)
        if nargin==0 || isempty(tpath) || ~ischar(tpath)
            tpath = uigetdir(opts.save_path,'Select folder for saving results:');
        end
        if tpath
            opts.save_path = tpath;
            h.text_save.Value = tpath;
        end
        figure(h.fig);
    end
    function selectAll(~,~)
        set([h.unreg,h.airway,h.scatnetAT,h.scatnetEmph,h.vessel,h.reg,h.prm,h.tprm],'Value',1);
    end
    function clearAll(~,~)
        set([h.unreg,h.airway,h.scatnetAT,h.scatnetEmph,h.vessel,h.reg,h.prm,h.tprm],'Value',0);
    end
    function preview(~,~)
        if isempty(preview_ind)
            warning('Select a row to preview.');
        else
            preview_path = C.DataPath{preview_ind};
            if isempty(cmiObj) || ~isvalid(cmiObj)
                cmiObj = CMIclass;
            end
            cmiObj.loadImg(0,preview_path);
        end
    end
    function selectCell(~,edata)
        preview_ind = edata.Indices(1) + find(page_ic,1) - 1;
    end
    function pipelineHELP(~,~)
        if isempty(h_help) || ~isvalid(h_help)
            ss = get(0,"ScreenSize");
            fn_help = fullfile(fileparts(which('cmi')),'scripts','LungCT','CTlung_Pipeline_Help.pdf');
            h_help = uifigure("Name","CTlung Pipeline Help","Position",[1,1,ss(3)/2,ss(4)]);
            uihtml(h_help,"Position",[0,0,h_help.Position(3:4)],"HTMLSource",fn_help);
        end
    end
    function saveCatalog(~,~)
        [fname,fpath] = uiputfile(opts.dcm_path,'Save Catalog');
        if fname
            fname = fullfile(fpath,fname);
            writetable(C,fname);
        end
    end
    function run(~,~)    
        % Check that each timepoint has Ins/Exp selected (and only one of each)
        if any(gp_valid == 3)
            warning('No case can have more than one Exp/Ins selected.');
        elseif any(gp_valid == 4)
            warning('Each case must have a unique UMlabel.');
        elseif ~isempty(C)
            % Run Pipeline on selected data
            
            if ~ismember('CaseNumber',C.Properties.VariableNames)
                C = addvars(C,ugroups_ic,'Before',1,'NewVariableNames',{'CaseNumber'});
            end
        
            % Set up output structure:
            empty_flag = false(ngroups,1);
            selected_data = struct('UMlabel',cell(1,ngroups),'StudyDate',cell(1,ngroups),'Scans',cell(1,ngroups));
            for ig = 1:ngroups
                g_ind = ugroups_ic==ig;
                tC = C( g_ind & ~cellfun(@isempty,C.Tag),:);
                if isempty(tC)
                    empty_flag(ig) = true;
                else
                    selected_data(ig).UMlabel = tC.UMlabel{1};
                    selected_data(ig).StudyDate = tC.StudyDate{1};
                    selected_data(ig).Scans = tC;
                end
            end
            selected_data(empty_flag) = [];

            CTlung_Pipeline_run(selected_data,opts);
            
        end


    end
end
