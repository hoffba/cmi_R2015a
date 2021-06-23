function [selected_data,home_path] = catalog_select(C)
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

% Put SeriesDescription column first:
[~,I] = sort(strcmp(colnames,'SeriesDescription'),'descend');
colnames = colnames(I);
C = C(:,I);

% Add Tag column:
[nscans,nfields] = size(C);
nfields = nfields+1;
colnames = [{'Tag'},colnames];
C = [repmat({''},nscans,1),C];

% Find unique Patient IDs:
ci = find(strcmp(colnames,'PatientName'),1);
[unames, ~, unames_ic] = unique(C(:,ci));
nnames = length(unames);
ntp = zeros(1,nnames);

for inames = 1:nnames
    selected_data(inames).PatientID = unames{inames};
    
    % Find unique timepoints:
    ci = find(strcmp(colnames,'StudyDate'),1);
    tC = C(unames_ic==inames,:);
    [utp, ~, utp_ic] = unique(cellfun(@num2str,tC(:,ci),'UniformOutput',false));
    ntp(inames) = length(utp);
    for itp = 1:ntp(inames)
        selected_data(inames).timepoint(itp).StudyDate = utp{itp};
        ttC = tC(utp_ic==itp,:);
        
        % Sort by SeriesNumber
        [~,I] = sort([ttC{:,strcmp(colnames,'SeriesNumber')}]);
        ttC = ttC(I,:);
        
        % Try to find Ins/Exp tags:
        ci = strcmp(colnames,'SeriesDescription');
        ie = find(contains(ttC(:,ci),'exp','IgnoreCase',true),1);
        if ~isempty(ie)
            ttC{ie,1} = 'Exp';
        end
        ie = find(contains(ttC(:,ci),'ins','IgnoreCase',true),1);
        if ~isempty(ie)
            ttC{ie,1} = 'Ins';
        end
        selected_data(inames).timepoint(itp).Scans = ttC;
        checkValid(inames,itp);
    end
end

%% Set up figure:
i_pid = 1;
i_tp = 1;
scsz = get(0,'screensize'); % screen size
fwidth = scsz(3)/2;
fheight = scsz(4)/3;
listwidth = 100; 
gap = 5;
pos = round([scsz(3)/4, scsz(4)/2, fwidth, fheight]);
hf = uifigure('Position',pos,'Name','Select DICOMcatalog data for processing:',...
    'CloseRequestFcn',@cancel_callback);
% guidata(hf,selected_data);
uilabel(hf,'Position',[ 1 , fheight-20 , listwidth , 20 ],...
    'Text','PatientID:','HorizontalAlignment','center');
uilabel(hf,'Position',[ listwidth+gap , fheight-20 , listwidth , 20 ],...
    'Text','StudyDate:','HorizontalAlignment','center');
htext_id = uitextarea(hf,'Position',[ 1 , fheight-40 , listwidth , 20 ],...
    'Value',unames{1},'ValueChangedFcn',@updatePatientID);
hlist_id = uicontrol(hf,'Style','list','Position',[ 1 , 1 , listwidth , fheight-40-gap ],...
    'Items',unames,'ValueChangedFcn',@selectPatient);
hlist_tp = uicontrol(hf,'Style','list','Position',[ listwidth+gap , 1 , listwidth , fheight-40-gap ],...
    'Items',{selected_data(i_pid).timepoint(:).StudyDate},'ValueChangedFcn',@selectTP);
colForm = repmat({''},1,nfields); colForm{1} = {' ','Exp','Ins'};
huit = uitable(hf,'Position',[ 2*(listwidth+gap) , 1, fwidth-2*(listwidth+gap) , fheight-20 ],...
    'ColumnName',colnames,'Data',selected_data(i_pid).timepoint(i_tp).Scans,...
    'ColumnFormat',colForm,'ColumnEditable',[true,false(1,nfields-1)],'CellEditCallback',@editTag);
uibutton(hf,'Position',[fwidth-50,fheight-20,50,20],...
    'Text','Done','BackgroundColor','green','ButtonPushedFcn',@done_callback);
uibutton(hf,'Position',[fwidth-100-gap,fheight-20,50,20],...
    'Text','Cancel','BackgroundColor','red','ButtonPushedFcn',@cancel_callback);

%% End script when window closes
waitfor(hf);

%% Set up output structure:
if ~isempty(selected_data)
    for iname = 1:nnames
        ntp = length(selected_data(iname).timepoint);
        for itp = 1:ntp
            
            % Fix empty tags (GUI requires non-empty selections, so it's shown as a space):
            ind = strcmp(selected_data(iname).timepoint(itp).Scans(:,1),' ');
            if any(ind)
                selected_data(iname).timepoint(itp).Scans{ind,1} = [];
            end
            
            % Convert Scans cell array to structure:
            selected_data(iname).timepoint(itp).Scans = ...
                cell2struct(selected_data(iname).timepoint(itp).Scans,colnames,2);
            
        end
    end
end

%% Set up callbacks:
    function selectPatient(src,~)
        i_pid = strcmp(src.Items,src.Value);
        % Update edit box:
        htext_id.Value = src.Value;
        % Update timepoints list:
        hlist_tp.Items = {selected_data(i_pid).timepoint(:).StudyDate};
        % Update scan table:
        i_tp = 1;
        huit.Data = selected_data(i_pid).timepoint(i_tp).Scans;
    end
    function updatePatientID(src,~)
        hlist_id.Items(strcmp(hlist_id.Items,hlist_id.Value)) = src.Value;
    end
    function selectTP(src,~)
        i_tp = strcmp(src.Items,src.Value);
        % Update scan table:
        huit.Data = selected_data(i_pid).timepoint(i_tp).Scans;
    end
    function editTag(~,event)
        selected_data(i_pid).timepoint(i_tp).Scans{event.Indices(1),event.Indices(2)} = event.NewData;
    end
    function checkValid(iname, itp)
        tags = {selected_data(iname).timepoint(itp).Scans(:,1)};
        n_exp = nnz(strcmp(tags,'Exp'));
        n_ins = nnz(strcmp(tags,'Ins'));
        selected_data(iname).timepoint(itp).valid = (n_exp==n_ins && ismember(n_exp,[0,1]));
    end
    function UDlistbox(iname, itp)
        str_id = hlist_id.Items;
        str_tp = hlist_tp.Items;
        if selected_data(iname).timepoint(itp).valid
        else
        end
    end
    function done_callback(~,~)
        % Check that each timepoint has Ins/Exp selected (and only one of each)
        errchk = false;
        for i = 1:nnames
            for j = 1:ntp(i)
                n_exp = nnz(strcmp(selected_data(i).timepoint(j).Scans(:,1),'Exp'));
                n_ins = nnz(strcmp(selected_data(i).timepoint(j).Scans(:,1),'Ins'));
                if ~(n_exp==n_ins && ismember(n_exp,[0,1]))
                    if ~errchk
                        warning('Each time point must have either nothing tagged or one of each Exp/Ins.');
                        warning('Invalid timepoints:');
                        errchk = true;
                    end
                    warning(['   ', selected_data(i).PatientID, ': ', ...
                        num2str(selected_data(i).timepoint(j).StudyDate)]);
                end
            end
        end
        if ~errchk
            hf.delete;
        end
    end
    function cancel_callback(~,~)
        selected_data = [];
        hf.delete;
    end
end
