% CMIcalss function
% Run custom script
function scriptRun(self,x,~)
spath = fileparts(which('cmi'));
if ~isempty(spath)
    str = struct2cell(dir(fullfile(spath,'scripts','*.m')));
    if ~isempty(str)
        [sel,ok] = listdlg('ListString',str(1,:));
        if ok
            if strcmp(get(x,'Tag'),'script_run')
                [~,str,~] = fileparts(str{1,sel});
                feval(eval(['@',str]),self);
            else % open file for editing
                edit(str{1,sel});
            end
        end
    else
        warndlg(['No script files found. (',spath,filesep,'scripts)'])
    end
end
