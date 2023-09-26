classdef ElxClass < handle
    properties (SetObservable, SetAccess=private, GetAccess=public)
        
        % Default for Linux
        elxdir = '';  % Directory for Elastix/Transformix
        xtstr = '';   % Directory for Xterm exe
        sepstr = '';  % Command input seperator
        
        % Stack of multi-parameter optimizations
        T0check = false; % Determines whether to apply initial transform during Elastix call
        Tx0guess = struct('i',{},'fname',{},'fpath',{}); % File names of initial guess transform (.txt)
        Tx0              % Initial transformation(s) (before optimization)
        Schedule = {};   % Schedule of coregistration parameters
        outfmt = '.nii' % Image output format
        sys = 0;
        
    end
    methods
        % Constructor
        function self = ElxClass
            if ismac
                fprintf('ElxClass : Mac\n');
                self.elxdir = '/usr/local/bin/';
                self.xtstr = '/opt/X11/bin/xterm -geometry 170x50 ';
                self.sepstr = ' ; ';
                self.sys = 1;
            elseif ispc
                fprintf('ElxClass : Windows\n');
                self.xtstr = 'cmd /c mode con: cols=100 lines=40 && ';
                self.sepstr = ' && ';
                self.sys = 2;
            else % Linux
                [stat,hostname] = system('uname -n');
                fprintf('ElxClass : Linux : hostname = %s\n',hostname);
                if contains(hostname,'.arc-ts.umich.edu')
                    % For Great Lakes, don't open a terminal
                    fprintf('Great Lakes!\n')
                    self.xtstr = 'module load Bioinformatics elastix ; module load gcc ; ';
                    self.sepstr = ' ; ';
                    self.sys = 4;
                elseif contains(hostname,'galban-ap-ps1a.med.umich.edu')
                    % Galban lab Tier2 server
                    fprintf('Galban Teir2!\n')
                    self.xtstr = '';
                    self.sepstr = ' ; ';
                    self.sys = 4;
                else
                    % Unknown system
                    fprintf('Unknown System!\n')
                    self.elxdir = '/opt/elastix/bin/';
                    self.xtstr = 'xterm -geometry 170x50 ';
                    self.sepstr = ' ; ';
                    self.sys = 3;
                end
            end
        end
    end
    methods (Static)
        pname = copyTfxChain(odir,fnames,itname)
    end
end