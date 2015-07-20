% RegClass function
function setTguess(self,x,edata)

if (nargin==1) || ischar(x) || (isa(x,'matlab.ui.control.UIControl') && strcmp(x.Tag,'button_Tguess'))
    if ischar(x)
        fname = x;
        fpath = '';
    else
        [fname,fpath] = uigetfile('*.txt','Select Transform Parameter File:');
    end
    if ischar(fname)
        fname = fullfile(fpath,fname);
        
        % Find all transforms in queue:
        tname = fname;
        allnames = {};
        while ~strcmp(tname,'NoInitialTransform') && exist(tname,'file')
            % Add to list of transforms:
            allnames = [allnames;{tname}];
            % Read file for next initial transform:
            fid = fopen(tname,'rt');
            tstr = fread(fid,'*char')';
            fclose(fid);
            tok = regexp(tstr,'\(InitialTransformParametersFileName \"(.*?)\"\)','tokens');
            tname = tok{1}{1};
        end
        
        % Update queue table:
        [fpath,fname,ext] = cellfun(@(x)fileparts(x),allnames,'UniformOutput',false);
        nt = length(fpath);
        C = [num2cell(true(nt,1)),strcat(fname,ext),fpath];
        set(self.h.table_Tguess,'Data',C);
        
        % Set ElxClass:
        self.elxObj.setT0guess(cell2struct(C',{'i','fname','fpath'}));
        
    end
elseif (isa(x,'matlab.ui.control.UIControl') && strcmp(x.Tag,'table_Tguess'))
    
    % Grab new selection value:
    i = edata.Indices(1);
    nval = edata.NewData;
    
    % Set new selection in ElxClass:
    str = self.elxObj.Tx0guess;
    str(i).i = nval;
    self.elxObj.setT0guess(str);
    
end

