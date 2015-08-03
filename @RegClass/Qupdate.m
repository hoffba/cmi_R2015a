% RegClass function
function Qupdate(self,x,edata)

% Read Qfile:
if exist(self.qfile,'file')
    fid = fopen(self.qfile,'rt');
    str = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    str = str{1};
else
    str = {};
end
nq = length(str);

% Extract names:
if nq
    if ispc
        expr = [' title (.*?)',self.sepstr];
    else
        expr = ' -T \"(.*?)\"';
    end
    enames = cellfun(@(x)regexp(x,expr,'tokens'),str);
    enames = [enames{:}]';
    enames(cellfun(@isempty,enames)) = {''};
else
    enames = {};
end

% Determine wheter a re-ordering was called from table:
if isa(x,'matlab.ui.control.Table') && nq
    if strcmp(x.Tag,'table_queue')
        
        % Grab new order value:
        i = edata.Indices(1);
        nval = round(str2double(edata.NewData));
        nstr = x.Data{i,2};
        
        % Double-check that Qfile line is still same:
        if ~strcmp(enames{i},nstr)
            warning('Queue file has changed, try again.')
        else
            
            % Re-order lines:
            if nval>nq
                nval = nq;
            elseif nval<0
                nval = 1;
            end
            tstr = str(i); tname = enames(i);
            str(i) = []; enames(i) = [];
            if nval==0
                % (0 means delete the line)
                nq = nq-1;
            else
                str = [str(1:nval-1);tstr;str(nval:end)];
                enames = [enames(1:nval-1);tname;enames(nval:end)];
            end
            
            % Overwrite the Qfile:
            fid = fopen(self.qfile,'wt');
            fprintf(fid,'%s\n',str{:});
            fclose(fid);
        end
    end
end

% Now update the displayed queue:
set(self.h.table_queue,'Data',...
    [ cellfun(@num2str,num2cell((1:nq)'),'UniformOutput',false), enames(:) , str(:) ]);
set(self.h.text_Qupdate,'String',datestr(now()));


