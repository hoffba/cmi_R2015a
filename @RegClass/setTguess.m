% RegClass function
function setTguess(self,x,edata)

if (nargin==1) || ischar(x) || (isa(x,'matlab.ui.control.UIControl') && strcmp(x.Tag,'button_TGload'))
    if ischar(x)
        [fpath,tname,ext] = fileparts(x);
        tname = [tname,ext];
    else
        [tname,fpath] = uigetfile('*.txt','Select Transform Parameter File:');
        if ~ischar(fpath)
            return;
        end
    end

    % Find all transforms in queue:
    allnames = {};
    tname = fullfile(fpath,tname);
    while ~strcmp(tname,'NoInitialTransform')
    % Add to list of transforms:
        allnames = [allnames;{tname}];
    % Read file for next initial transform:
        fid = fopen(tname,'rt');
        tstr = fread(fid,'*char')';
        fclose(fid);
    % Find file name
        tok = regexp(tstr,'\(InitialTransformParametersFileName \"(.*?)\"\)','tokens');
    % Assume all initial transform files are in same folder (corrects for copied files)
        [~,tname,ext] = fileparts(tok{1}{1});
        if ~isempty(ext)
            tname = fullfile(fpath,[tname,ext]);
        end
    end

    % Check that files exist:
    chk = cellfun(@(x)exist(x,'file'),allnames);
    if all(chk)
        % Update queue table:
        [fpath,fname,ext] = cellfun(@(x)fileparts(x),allnames,'UniformOutput',false);
        nt = length(fpath);
        C = [num2cell(true(nt,1)),strcat(fname,ext),fpath];
        set(self.h.table_Tguess,'Data',C);

        % Set ElxClass:
        self.elxObj.setT0guess(cell2struct(C',{'i','fname','fpath'}));
    else
        error([sprintf('Files were not found:\n'),sprintf('    %s\n',allnames{~chk})]);
    end
        
elseif isa(x,'matlab.ui.control.UIControl') && strcmp(x.Tag,'button_TGclear')
    self.elxObj.setT0guess(struct('i',{},'fname',{},'fpath',{}));
    set(self.h.table_Tguess,'Data',{});
elseif isa(x,'matlab.ui.control.Table') && strcmp(x.Tag,'table_Tguess')
    
    % Grab new selection value:
    i = edata.Indices(1);
    nval = edata.NewData;
    
    % Set new selection in ElxClass:
    str = self.elxObj.Tx0guess;
    str(i).i = nval;
    self.elxObj.setT0guess(str);
    
end

