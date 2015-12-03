% Function to generate listmode trigger times in Bruker .lst files
%       if forgot to start trigger and ".._pm.lst" was not created
function genPMlist(aqname,toff,dt)
% Inputs: fname = ".._aq.lst" full filename
%         toff  = offset time (seconds)
%         dt    = time between triggers

p = 1; % seconds

% GUI if called without inputs
if nargin==0
    [aqname,fpath] = uigetfile('*_aq.lst','Select AQ list file:');
    if ischar(aqname)
        aqname = fullfile(fpath,aqname);
        val = inputdlg({'Shift time by how much (seconds)?',...
                        'Time between triggers:'},...
                       'LST Options',1,{'0','1'});
        if isempty(val)
            reture;
        else
            toff = str2double(val{1});
            dt = str2double(val{2});
        end
    else
        return;
    end
elseif ~exist(aqname,'file') || ~strcmp(aqname(end-6:end),'_aq.lst')
    error('Invalid input file');
else
    if nargin<2
        toff = 0;
    end
    if nargin<3
        dt = 1;
    end
end
    
% Set output name
pmname = aqname;
pmname(end-5:end-4) = 'pm';

% Determine acquisition time limits:
t = readLST(aqname);
t0 = t(1);
t1 = t(end);

% Next, generate PM timestamps within acquisition limits:
t = ((t0+toff) : dt : (t1+dt));

% Save PM list:
saveLST(t,pmname);

