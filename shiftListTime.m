% Function to shift listmode times in Bruker .lst files
function shiftListTime(fname,val)

% GUI if called without inputs
if nargin==0
    [fname,path] = uigetfile('*pm.lst','Select PM list file:');
    if fname~=0
        fname = fullfile(path,fname);
        val = str2double(inputdlg('Shift time by how much (seconds)?',...
                                  'LST Time Shift',1,{'0'}));
    else
        fname = [];
    end
elseif (nargin~=2)
    val = [];
end

if ~(isempty(fname) || isempty(val)) && ischar(fname) && ...
        exist(fname,'file') && strcmp(fname(end-3:end),'.lst') && ...
        isnumeric(val) && ~isnan(val)
    
    % First, determine acquisition timestamps:
    [path,bname,~] = fileparts(fname);
    acqname = fullfile(path,[bname(1:end-2),'aq.lst']);
    aqtimes = [];
    if exist(acqname,'file')
        fid = fopen(acqname,'r');
        if fid>2
            str = strsplit(fread(fid,inf,'*char')',{'=','\r'});
            fclose(fid);
            nts = (length(str)-1)/2;
            aqtimes = zeros(1,nts);
            for i = 1:nts
                t = str2double(strsplit(str{1+2*i},':'));
                aqtimes(i) = (t(1)*60 + t(2))*60 + t(3) + t(4)/1000;
            end
        end
    end
    
    % Next, adjust PM timestamps:
    if ~isempty(aqtimes)
        pmtimes = [];
        tname = fullfile(path,[bname,'_orig.lst']);
        if ~exist(tname,'file')
            copyfile(fname,tname);
        end
        
        % Read in PM timestamps:
        fid = fopen(fname,'r');
        if fid>2
            str = fread(fid,inf,'*char');
            fclose(fid);
            eqind = strfind(str,'=');
            strspl = strsplit(str',{'=','\r'});
            nts = length(eqind)-1;
            pmtimes = zeros(1,nts);
            for i = 1:nts
                t = str2double(strsplit(strspl{3+2*i},':'));
                pmtimes(i) = (t(1)*60 + t(2))*60 + t(3) + t(4)/1000 + val;
            end
        end
        
        % Write adjusted timestamps to file:
        if ~isempty(pmtimes)
            
            % Adjust timestamps:
            % Need to make sure that there's one PM timestamp both before
            %       the first and after the last AQ timestamp.
            avgdt = mean(diff(pmtimes));
            n = nnz(pmtimes < aqtimes(1));
            if n>1 % need to remove pre-img times
                pmtimes(1:n-1) = [];
            elseif n==0 % need to add pre-img times
                n = ceil((pmtimes(1) - aqtimes(1))/avgdt);
                pmtimes = [pmtimes(1)-(n:1)*avgdt , pmtimes];
            end
            n = nnz(pmtimes >= aqtimes(end));
            if n>1 % need to remove post-img times
                pmtimes(end-n:end) = [];
            elseif n==0 % need to add post-img times
                n = ceil((pmtimes(end) - aqtimes(end))/avgdt);
                pmtimes = [pmtimes , pmtimes(end)+(1:n)*avgdt];
            end
            
                
            
            % Write adjusted pm.lst file:
            fid = fopen(fullfile(path,[bname,'_MOD.lst']),'w');
            if fid>2
                if strcmp(str{1},'[PM timemarks]') && strncmp(str{2},'total=',6)
                    npm = str2double(strtok(str{2},'total='));
                    fprintf(outfid,'%s\n',str{1},str{2});
                    hw = waitbar(0);
                    for i = 1:npm
                        str = fgets(fid);
                        strs = strsplit(str,'=');
                        t = str2double(strsplit(strs{2},':'));
                        tval = (t(1)*60 + t(2))*60 + t(3) + t(4)/1000 + val;
                        t(4) = round(mod(tval,1) * 1000);
                        tval = floor(tval);
                        t(3) = mod(tval,60);
                        tval = floor(tval/60);
                        t(2) = mod(tval,60);
                        t(1) = floor(tval/60);
                        str(1:17) = sprintf('%04u=%02u:%02u:%02u:%03u',i-1,t);
                        fprintf(fid,'%s',str);
                        waitbar(i/npm,hw,[num2str(i),' of ',num2str(npm),' done'])
                    end
                    delete(hw);
                    fclose(fid);
                    fclose(fid);
                else
                    error('Invalid file.');
                end
            else
                error('File could not be opened.');
            end
        end
    end
end

