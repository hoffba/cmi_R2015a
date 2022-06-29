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
opts = struct('dcmpath','',...
              'save_path','',...
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
C = [];

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

%% Set catalog display:
opts.dcmpath = selectCatalog(opts.dcmpath);
figure(h.fig);

%% Set path for saving results:
setSavePath(opts.save_path);

%% End script when window closes
waitfor(h.fig);

%% Gather outputs
if ~isempty(C)
    %% Determine tags:
    tagstr = repmat({''},size(C,1),1);
    tagstr(C.Exp) = {'Exp'};
    tagstr(C.Ins) = {'Ins'};
    C = addvars(C,ugroups_ic,tagstr,'Before',1,'NewVariableNames',{'CaseNumber','Tag'});
    
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
    fname = opts.dcmpath;
    if ~isempty(fname)
        C = removevars(C,{'Ins','Exp'});
        if isfolder(fname)
            fname = fullfile(fname,'DICOMcatalog_select.csv');
        end
        writetable(C,fname);
    end

end

%% Set up callbacks:
    function validateInputs(inputs)
        p = inputParser;
        addParameter(p,'dcm_path','',@(x)ischar(x)&&isfolder(x));
        addParameter(p,'save_path','',@(x)ischar(x)&&isfolder(x));
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
            nfields = length(colnames);
            
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
            if ismember('CaseNumber',colnames)
                ugroups_ic = C.CaseNumber;
            else
                C(cellfun(@isempty,C.UMlabel),:) = [];
                C = sortrows(C,{'StudyDate','StudyID','SeriesNumber'});

                % Find groups of scans with unique PatientName and StudyDate
                if isnumeric(C.StudyDate)
                    C.StudyDate = cellfun(@num2str,num2cell(C.StudyDate),'UniformOutput',false);
                end
                [~,~,ugroups_ic] = unique(strcat(C.PatientName,C.StudyDate,C.StudyID));
            end
            ngroups = numel(unique(ugroups_ic));

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
            checkValid;
            
        end
        h.text_cat.Value = str;
    end
    function setSavePath(tpath,~)
        if nargin==0 || isempty(tpath) || ~ischar(tpath)
            tpath = uigetdir('Select folder for saving results:');
        end
        if ~isempty(tpath)
            opts.save_path = tpath;
            h.text_save.Value = tpath;
        end
    end
    function selectAll(~,~)
        set([h.unreg,h.airway,h.scatnet,h.vessel,h.reg,h.prm,h.tprm],'Value',1);
    end
    function clearAll(~,~)
        set([h.unreg,h.airway,h.scatnet,h.vessel,h.reg,h.prm,h.tprm],'Value',0);
    end
end
