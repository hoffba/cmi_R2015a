function [selected_data,fpath] = catalog_select_2(varargin)
% GUI for selection of DICOM data from DICOMcatalog information
%   For use in Ins/Exp CT analysis, looking for two matching scans per time
% Input: C = table of DICOM information
%          OR = filename of catalog .csv to read in
% Output: selected_data = structure containig DICOM info for selected data
%                      .PatientID = patient name
%                      .timepoint = structure array containing timepoint info
%                                .StudyDate = string date
%                                .Scans = structure array with scan info

C = [];
fpath = '';
fname = '';
cat_flag = false; % T/F run dcmCatalog

if nargin==0
    fpath = uigetdir('Select folder for processing.');
    if ~fpath
        return;
    end

    % Find existing DICOMcatalog files:
    fnames = dir(fullfile(fpath,'*DICOMcatalog*.csv'));
    answer = listdlg('PromptString','Select a catalog:',...
        'ListString',{'* Create New *',fnames.name});
    if isempty(answer)
        return;
    elseif answer == 1
        cat_flag = true;
    elseif answer > 1
        fname = fullfile(fnames(answer-1).folder,fnames(answer-1).name);
        fpath = '';
    end
else
    for ii = 1:nargin
        val = varargin{ii};
        if istable(val)
            C = val;
            cat_flag = false;
            break;
        elseif ischar(val) && isfolder(val)
            fpath = val;
            cat_flag = true;
        elseif ischar(val) && strcmp(val(end-3:end),'.csv')
            fname = val;
            fpath = fileparts(fname);
        else
            warning('Invalid input');
            return;
        end
    end
end

if cat_flag || ~exist(fname,'file')
    C = dcmCatalog(fpath,fname);
elseif ~isempty(fname) && exist(fname,'file')
    iopt = detectImportOptions(fname);
    C = readtable(fname,iopt);
else
    error('Invalid input string. Must be path to DICOMcatalog.csv file.');
end

%% Validate table input:
colnames = C.Properties.VariableNames;
req_fields = {'SeriesDescription','PatientName','StudyDate','SeriesNumber',...
              'ConvolutionKernel','SliceThickness','Slices'};
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
nscans = size(C,1);
C = sortrows(C,{'UMlabel','StudyDate','StudyID','SeriesNumber'});

% Find groups of scans with unique PatientName and StudyDate
if isnumeric(C.StudyDate)
    C.StudyDate = cellfun(@num2str,num2cell(C.StudyDate),'UniformOutput',false);
end
[~,~,ugroups_ic] = unique(strcat(C.PatientName,C.StudyDate,C.StudyID));
ngroups = max(ugroups_ic);


%% Set up figure:
scsz = get(0,'screensize'); % screen size
fwidth = scsz(3)/2;
fheight = scsz(4)*2/3;
gap = 5;
hf = uifigure('Position',round([(scsz(3)-fwidth)/2, scsz(4)/6, fwidth, fheight]),...
    'Name','Select DICOMcatalog data for processing:');%,'CloseRequestFcn',@cancel_callback);

colEdit = false(1,nfields); colEdit([1,2,3]) = true;
colWidth = repmat({'auto'},1,nfields); colWidth(1:2) = {40,40};
huit = uitable(hf,'Position',[ 1 , 20 , fwidth , fheight-110 ], 'Data',C,...
    'ColumnEditable',colEdit,'CellEditCallback',@editCell,'ColumnWidth',colWidth);

% dtypes = cell(1,nfields-2);
% for i_fld = 1:(nfields-2)
%     dtypes{i_fld} = class(C.(colnames{i_fld+2}));
% end
F = table('Size',[2,nfields-2],'VariableTypes',repmat({'cellstr'},1,nfields-2),'VariableNames',[{'Tag'},colnames(4:end)]);
F.Tag = {'Exp';'Ins'};
htfilt = uitable(hf,'Position',[ 1 , fheight-90 , fwidth , 90 ],...
    'Data',F,'ColumnEditable',[false,true(1,nfields-2)],'CellEditCallback',@applyFilter);

uibutton(hf,'Position',[1,1,50,20],'Text','Clear','BackgroundColor','blue','ButtonPushedFcn',@clearTags);
uibutton(hf,'Position',[fwidth-50,1,50,20],...
    'Text','Done','BackgroundColor','green','ButtonPushedFcn',@done_callback);
uibutton(hf,'Position',[fwidth-100-gap,1,50,20],...
    'Text','Cancel','BackgroundColor','red','ButtonPushedFcn',@cancel_callback);

gp_valid = zeros(ngroups,1);

pause(1);
applyFilter;

%% End script when window closes
waitfor(hf);

%% Set up output structure:
if isempty(C)
    selected_data = [];
else
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
    function applyFilter(~,eventdata)
        
        if nargin==0 % filter all
            rowi = 1:2;
        elseif nargin==2 % filter single
            rowi = eventdata.Indices(1);
        end
        
        % Set tags:
        for i = rowi
            if any(~cellfun(@isempty,table2cell(htfilt.Data(i,2:end))))
                TF = true(nscans,1);
                for j = 2:numel(htfilt.ColumnName)
                    tname = htfilt.ColumnName{j};
                    if iscell(htfilt.Data.(tname)) && ~isempty(htfilt.Data.(tname){i})
                        TF = TF & contains(C.(tname),htfilt.Data.(tname){i});
                    elseif isnumeric(htfilt.Data.(tname)) && htfilt.Data.(tname)(i)
                        TF = TF & (C.(tname)==htfilt.Data.(tname)(i));
                    end
                end
            else
                TF = false(nscans,1);
            end
            huit.Data.(htfilt.Data{i,1}{1}) = TF;
        end

        checkValid;
    end
    function editCell(hObject,eventdata)
        if ismember(eventdata.Indices(2),[1,2])
            % Edited Tag
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
            n_exp = nnz(huit.Data.Exp(irow));
            n_ins = nnz(huit.Data.Ins(irow));
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
            addStyle(huit,uistyle('BackgroundColor',tcolor),'row',irow);
        end
    end
    function clearTags(~,~)
        huit.Data(:,1:2) = {false};
        checkValid;
    end
    function done_callback(hObject,~)
        % Check that each timepoint has Ins/Exp selected (and only one of each)
        if any(gp_valid>2)
            warning('Each time point must have either nothing tagged or one of each Exp/Ins.');
        else
            C = huit.Data;
            hObject.Parent.delete;
        end
    end
    function cancel_callback(hObject,~)
        C = [];
        hObject.Parent.delete;
    end
end
