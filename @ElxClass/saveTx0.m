% ElxClass function
% Save Elastix Initial Transform file
function fname = saveTx0(self,x,~)

fname = '';
if ~isempty(self.Tx0)
    if (nargin==3) && ishandle(x)
        [fname,path] = uiputfile('*.txt','Save Initial Transform',...
                                         'InitialTransform.txt');
        if fname~=0
            fname = fullfile(path,fname);
        else
            fname = '';
        end
    elseif ~ischar(x)
        fname = '';
    else
        fname = x;
    end
    if ~isempty(fname)
        fid = fopen(fname,'w');
        if fid>2
            fstr = fieldnames(self.Tx0);
            for i = 1:length(fstr)
                val = self.Tx0.(fstr{i});
                vstr = '';
                if ischar(val)
                    vstr = [' "',val,'"'];
                elseif iscellstr(val)
                    vstr = sprintf(' "%s"',val{:});
                elseif all(round(val)==val)
                    vstr = sprintf(' % .0f',val);
                elseif isnumeric(val)
                    vstr = sprintf(' % .8f',val);
                end
                fprintf(fid,'(%s%s)\n',fstr{i},vstr);
            end
            fclose(fid);
        else
            error(['Could not open file: ',fname])
        end
    end
end