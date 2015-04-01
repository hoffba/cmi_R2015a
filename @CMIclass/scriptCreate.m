% CMIclass function
function scriptCreate(~,~,~)
str = inputdlg('Input script name:','Create Script',1);
spath = fullfile(fileparts(which('cmi')),'scripts');
if ~(isempty(str) || isempty(spath))
    fname = fullfile(spath,[str{1},'.m']);
    if exist(fname,'file')
        error('File already exists!')
    else
        % Initialize m-file
        fid = fopen(fname,'w');
        fprintf(fid,'%s\n',...
                '% CMI script',...
               ['function ',str{1},'(cmiObj)'],...
                '%   Input: cmiObj = CMIclass object containing current settings',...
                '');
        fclose(fid);
        % Open file for editing
        edit(fname);
    end
end