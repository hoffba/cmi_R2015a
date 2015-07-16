function batch_convertImages(ext)
% Batch function for converting images to different types

if nargin==0
    ext = '.mhd';
end

% Select files for conversion
fnames = {};
go = true;
while go
    [fname,fpath] = uigetfile('*','Select image for conversion:','MultiSelect','on');
    if iscellstr(fname)
        fnames = [fnames;cellfun(@(x)fullfile(fpath,x),fname,'UniformOutput',false)'];
    elseif ischar(fname)
        fnames = [fnames;{fullfile(fpath,fname)}];
    else
        go = false;
    end
end

if ~isempty(fnames)

    % GUI for renaming files:
    hf = figure('Units','normalized','ToolBar','none','MenuBar','none');
    ht = uitable(hf,'Data',[fnames,fnames],...
                    'Units','normalized',...
                    'Position',[0.02,0.1,0.96,0.88],...
                    'ColumnName',{'From:(.img)','To:(.mhd)'},...
                    'ColumnEditable',[false,true],...
                    'ColumnWidth',{200,200});
    hb = uicontrol(hf,'Style','pushbutton',...
                      'Units','normalized',...
                      'Position',[0.5,0.02,0.48,0.08],...
                      'String','OK',...
                      'Callback',@(hObject,~)delete(hObject));
    waitfor(hb); % Continue when OK button is clicked
    
    % Grab names from uitable:
    nnames = ht.Data(:,2);
    close(hf);

    % Create/Select output directory:
    odir = uigetdir(pwd,'Select directory for saving converted files:');

    % Loop over files:
    img = ImageClass;
    hw = waitbar(0);
    for i = 1:nf
        waitbar(i/nf,hw,['Converting image: ',nnames{i}])
        if exist(fnames{i},'file')
            img.loadImg(0,fnames(i));
            img.saveImg(1,fullfile(odir,nnames{i}),1);
        else
            warning(['File does not exist: ',fnames{i}]);
        end
    end
    delete(hw);
end