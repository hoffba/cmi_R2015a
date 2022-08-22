function stat = check_data(self,varargin)
% Returns requested data matrices

flds = fieldnames(self.dat);
for i = 1:numel(varargin)
    tag = varargin{i};
    if ismember(tag,flds)
        if isempty(self.dat.(tag).mat)
            
            % Determine filename
            fname = fullfile(self.procdir,[self.fn_base,'.',
            
            if exist(fname,'file')
                
                % Read from file:
                self.writelog('Reading from file: %s\n',fname);
                info = niftiinfo(fname);
                
                % Check image info match
                stat = self.check_info(info,tag);
                
                if stat
                    self.dat.(tag).mat = niftiread(fname);
                end
                
            else
                
                % Need to run process to generate requested data:
                ind = find(cellfun(@(x)ismember(tag,x),self.procs(:,2)),1);
                if isempty(ind)
                    warning('Requested data not defined: %s',tag');
                else
                    self.run(self.procs{ind,1});
                end
                
            end
        end
        dat_out = self.dat.(tag).mat;
    else
        warning('Invalid input tag: %s',tag);
    end
    varargout{i} = dat_out;
end