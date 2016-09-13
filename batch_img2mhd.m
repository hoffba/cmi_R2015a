function batch_img2mhd
% Batch function for converting/renaming Varian MRI data to MHD format.
% Creates new directory in selected study directory named MHD for saving
%   images.

% Find all .img folders in a user-selected folder
fdir = uigetdir(pwd);
fnames = dir(fullfile(fdir,'*.img'));
fnames = {fnames(:).name}';

if isempty(fnames)
    warning('No .img folders found in selected directory.')
else
    % Don't convert scout images
    fnames(strncmp('scout',fnames,5)) = [];
    nf = length(fnames);
    
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
    
    % Create saving directory:
    odir = fullfile(fdir,'MHD');
    if ~exist(odir,'dir')
        stat = mkdir(odir);
        if ~stat
            error('Could not create ./MHD directory. Check permissions.')
        end
    end
    
    % Loop over files:
    img = ImageClass;
    hw = waitbar(0);
    for i = 1:nf
        waitbar(i/nf,hw,['Converting image: ',nnames{i}])
        imgdir = fullfile(fdir,fnames{i});
        tfname = dir(fullfile(imgdir,'*.fdf'));
        if isempty(tfname)
            warning(['No .fdf files found in .img directory: ',imgdir])
        else
            img.loadImg(0,{fullfile(imgdir,tfname(1).name)});
            if (img.dims(4)==4) && ~isempty(strfind(imgdir,'_iadc_'))
                tvec = [2,4];
            else
                tvec = 1:img.dims(4);
            end
            [~,tfname] = fileparts(nnames{i});
            img.saveImg(tvec,fullfile(odir,[tfname,'.mhd']),1);
        end
    end
    delete(hw);
end

