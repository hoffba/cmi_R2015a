function [selected_data,home_path] = catalog_select_2(C)
% GUI for selection of DICOM data from DICOMcatalog information
%   For use in Ins/Exp CT analysis, looking for two matching scans per time
% Input: C = cell array of DICOM information
% Output: selected_data = structure containig DICOM info for selected data
%                      .PatientID = patient name
%                      .timepoint = structure array containing timepoint info
%                                .StudyDate = string date
%                                .Scans = structure array with scan info

if nargin==0 || isempty(C)
    % Use GUI to select DICOMcatalog.csv
    [tfname,home_path] = uigetfile('*.csv','Select DICOMcatalog:');
    if ~home_path
        selected_data = [];
        return
    end
    C = fullfile(home_path,tfname);
end
if ischar(C)
    if exist(C,'file') && strcmp(C(end-3:end),'.csv')
        C = cmi_csvread(C);
    else
        error('Invalid input string. Must be path to DICOMcatalog.csv file.');
    end
elseif ~iscell(C)
    error('Invalid input type.');
else
    home_path = pwd;
end

%% Validate cell array input:
% Separate header from scan info:
colnames = C(1,:);
C(1,:) = [];
req_fields = {'SeriesDescription','PatientName','StudyDate','SeriesNumber'};
i_member = ismember(req_fields,colnames);
if ~all(ismember({'SeriesDescription','PatientName','StudyDate','SeriesNumber'},colnames))
    error(['Missing required fields: ',strjoin(req_fields(~i_member),', ')]);
end

%% Set up display data:

% Re-order columns for easier reading:
I = cellfun(@(x)find(strcmp(colnames,x),1),req_fields);
tI = 1:length(colnames);
tI(I) = [];
I = [I tI];
colnames = colnames(I);
C = C(:,I);

% Add Tag column:
[nscans,nfields] = size(C);
nfields = nfields+1;
colnames = [{'Tag'},colnames];
C = [repmat({''},nscans,1),C];

% Sort array
ci = [ find(strcmp(colnames,'PatientName')) ,...
       find(strcmp(colnames,'StudyDate')) ,...
       find(strcmp(colnames,'SeriesNumber')) ];
C = sortrows(C,ci);

% Find groups of scans with unique PatientName and StudyDate
ci_PatientName = find(strcmp(colnames,'PatientName'),1);
ci_SeriesDescription = find(strcmp(colnames,'SeriesDescription'),1);
ci_StudyDate = find(strcmp(colnames,'StudyDate'),1);
[~,~,unames_ic] = unique(C(:,ci_PatientName));
[~,~,ugroups_ic] = unique([unames_ic,[C{:,ci_StudyDate}]'],'rows');
ngroups = max(ugroups_ic);

% Loop over scan groups
gp_valid = true(ngroups,1);
for i_gp = 1:ngroups

    % Try to find Ins/Exp tags:
    idx = find(ugroups_ic==i_gp);
    tC = C(idx,ci_SeriesDescription);
    gp_valid(i_gp) = true;
    i_exp = find(contains(tC,'exp','IgnoreCase',true),1);
    if ~isempty(i_exp)
        C{idx(i_exp),1} = 'Exp';
        gp_valid(i_gp) = ~gp_valid(i_gp);
    end
    i_ins = find(contains(tC,'ins','IgnoreCase',true),1);
    if ~isempty(i_ins)
        C{idx(i_ins),1} = 'Ins';
        gp_valid(i_gp) = ~gp_valid(i_gp);
    end
        
end

%% Set up figure:
scsz = get(0,'screensize'); % screen size
fwidth = scsz(3)/2;
fheight = scsz(4)*2/3;
gap = 5;
hf = uifigure('Position',round([(scsz(3)-fwidth)/2, scsz(4)/6, fwidth, fheight]),...
    'Name','Select DICOMcatalog data for processing:','CloseRequestFcn',@cancel_callback);
colForm = repmat({''},1,nfields); colForm{1} = {' ','Exp','Ins'};
colEdit = false(1,nfields); colEdit([1,3]) = true;
huit = uitable(hf,'Position',[ 1 , 20 , fwidth , fheight-20 ],'ColumnName',colnames,...
    'Data',C,'ColumnFormat',colForm,'ColumnEditable',colEdit,'CellEditCallback',@editCell);
uibutton(hf,'Position',[fwidth-50,1,50,20],...
    'Text','Done','BackgroundColor','green','ButtonPushedFcn',@done_callback);
uibutton(hf,'Position',[fwidth-100-gap,1,50,20],...
    'Text','Cancel','BackgroundColor','red','ButtonPushedFcn',@cancel_callback);

bgs_white = uistyle('BackgroundColor','white');
bgs_gray = uistyle('BackgroundColor','#cecece');
bgs_red = uistyle('BackgroundColor','#ffc4c4');

checkValid(huit);

%% End script when window closes
waitfor(hf);

%% Set up output structure:
if isempty(C)
    selected_data = [];
else
    selected_data = struct('PatientName',cell(1,ngroups),'StudyDate',cell(1,ngroups),'Scans',cell(1,ngroups));
    for ig = 1:ngroups
        g_ind = ugroups_ic==ig;
        tC = C( g_ind & ismember(C(:,1),{'Exp','Ins'}) ,:);
        if ~isempty(tC)
            selected_data(ig).PatientName = C{find(g_ind,1),3};
            selected_data(ig).StudyDate = num2str(C{find(g_ind,1),4});
            selected_data(ig).Scans = cell2struct(tC, colnames , 2);
        end
    end
end

%% Set up callbacks:
    function editCell(hObject,eventdata)
        if eventdata.Indices(2)==1
            % Edited Tag
            checkValid(hObject);
        else 
            % Edited PatientID
            if isempty(regexp(eventdata.NewData , '[/\*:?"<>|]', 'once'))
                % Set ID for all scans in this case group:
                hObject.Data(ugroups_ic==ugroups_ic(eventdata.Indices(1)),3) = {eventdata.NewData};
            else
                warning('Invalid PatientID, try again.');
                hObject.Data{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.PreviousData;
            end
        end
    end
    function checkValid(hObject,~)
        for i = 1:ngroups
            irow = find(ugroups_ic==i);
            tags = hObject.Data(irow,1);
            n_exp = nnz(strcmp(tags,'Exp'));
            n_ins = nnz(strcmp(tags,'Ins'));
            gp_valid(i) = (n_exp==n_ins && ismember(n_exp,[0,1]));
            if ~gp_valid(i)
                addStyle(hObject,bgs_red,'row',irow);
            elseif rem(i,2)
                addStyle(hObject,bgs_white,'row',irow);
            else
                addStyle(hObject,bgs_gray,'row',irow);
            end
        end
    end
    function done_callback(hObject,~)
        % Check that each timepoint has Ins/Exp selected (and only one of each)
        if any(~gp_valid)
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