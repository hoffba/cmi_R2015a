% Function to generate listmode trigger times in Bruker .lst files
%       if forgot to start trigger and ".._pm.lst" was not created
function genPMlist(fname,val)
% Inputs: fname = ".._aq.lst" full filename
%         val   = offset time (seconds)

p = 1; % seconds

% GUI if called without inputs
if nargin==0
    [fname,path] = uigetfile('*_aq.lst','Select AQ list file:');
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
        exist(fname,'file') && strcmp(fname(end-6:end),'_aq.lst') && ...
        isnumeric(val) && ~isnan(val)
    
    % Determine file names
    aqname = fname;
    [path,fname,~] = fileparts(fname);
    pmname = fullfile(path,[fname(1:end-2),'pm.lst']);
    
    % First, determine acquisition timestamps:
    if exist(aqname,'file')
        fid = fopen(aqname,'r');
        if fid>2
            str = strsplit(fread(fid,inf,'*char')',{'=','\r'});
            fclose(fid);
            t0 = str2t(str{3});
            t1 = str2t(str{floor((length(str)-1)/2)*2 + 1});
        end
    end
    
    % Modulate time offset by the trigger period
    val = mod(val,p);
    
    % Next, generate PM timestamps:
    pmtimes = ((t0-4*p) : p : (t1+4*p)) + val;
    npm = length(pmtimes);
    fid = fopen(pmname,'w');
    if fid>2
        fprintf(fid,'[PM timemarks]\rtotal=%04u\r',npm);
        for i = 1:npm
            fprintf(fid,'%04u=%s\r',i-1,t2str(pmtimes(i)));
        end
        fclose(fid);
    end
end

function t = str2t(str)
% Time string: HH:MM:SS:ddd
% Returns time in seconds
tsep = str2double(strsplit(str,':'));
t = (tsep(1)*60 + tsep(2))*60 + tsep(3) + tsep(4)/1000;

function str = t2str(tval)
t(4) = round(mod(tval,1) * 1000);
tval = floor(tval);
t(3) = mod(tval,60);
tval = floor(tval/60);
t(2) = mod(tval,60);
t(1) = floor(tval/60);
str = sprintf('%02u:%02u:%02u:%03u',t);
