% ElxClass function
% Save Elastix Parameter file
function fname = savePar(self,ind,fname)
% Save Elastix Parameter File (.txt)
% fname = savePar();
%           Calls GUI to select file name/directory
% fname = savePar([ind],'path/fname.txt');
% fname = savePar([ind],{'path/fname.txt','path/fname.txt',...});

ns = length(self.Schedule);
if ns~=0
    
    % Parse inputs for Schedule indices and output file names
    if (nargin==1)
        [fname,path] = uiputfile('*.txt','Save Elastix Parameters',...
                                         'ElastixParameters.txt');
        if fname~=0
            ind = 1:ns;
            fname = fullfile(path,fname);
        else
            fname = {};
        end
    elseif (nargin<3) || ~(iscellstr(fname)||ischar(fname))
        fname = {};
    end
    if ~isnumeric(ind) || isempty(ind)
        error('First input must be numerical indices within length(Schedule)')
    end
    
    if ~isempty(fname)
        ni = length(ind);
        if ischar(fname)
            if ni==1
                fname = {fname};
            else
                fname = strcat(fname(1:end-4),'_',...
                               cellfun(@num2str,num2cell(ind),'UniformOutput',false),...
                               '.txt');
            end
        end
        if length(fname)<ni
            error('Number of input file names does not match indices.');
        end
        
        % Loop over all Schedule steps:
        for istep = 1:ni
            fid = fopen(fname{istep},'w');
            if fid>2
                fstr = fieldnames(self.Schedule{ind(istep)});
                for i = 1:length(fstr)
                    val = self.Schedule{ind(istep)}.(fstr{i});
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
                error(['Could not open file: ',tname])
            end
        end
    end
end