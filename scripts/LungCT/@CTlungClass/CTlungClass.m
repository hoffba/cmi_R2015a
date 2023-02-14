classdef CTlungClass < handle
    properties (SetObservable, SetAccess=private, GetAccess=public)
        
        % Basic case info
        gui_flag = false;   % Switch for using GUI prompts
        id = "";        % Subject ID
        procdir = "";   % Processing directory
        elxdir = "";    % Elastix directory in procdir
        fn_log = "";    % File name of processing log
        dat % Data structure containing relevant pipeline information:
            % i.e. filename tags, loaded image data, info, options
            % (See method: initializeDat)
        
        % Pipeline options
        ext =           '.nii.gz';  % format for saving data
        gl_flag =       true;       % check to run on Great Lakes
        username =      '';         % user's UMICH uniquename
        orient_flag =   true;       % auto-check for image orientation
        swap_flag =     true;       % Flag for EXP/INS based on volume
        seg_method =    4;          % Segmentation method - 4=YACTA
        unreg_flag =    true;       % Flag to run unreg analysis
        airway_flag =   true;       % Flag to run airway analysis (YACTA)
        scatnet_flag =  true;       % Flag to perform scatnet (LAA)
        vessel_flag =   true;       % Flag to run vessel analysis
        reg_flag =      true;       % Flag to perform registration
        quickreg =      false;      % Skip last resolution of warp
        prm_flag =      true;       % Flag to run PRM analysis
        tprm_flag =     true;       % Flag to run tPRM analysis
        
        % Analysis settings
        yacta = '';                     % YACTA options
        prm = struct('check',true,...   % PRM settings
                     'thresh',[],...
                     'cutoff',[],...
                     'cmap',[],...
                     'prmmap',[],...
                     'SPopts',[]);
        tprm = struct('check',true,...  % tPRM settings
                      'thresh',1:4,...
                      'tmode','==',...
                      'n',10*ones(1,3),...
                      'gridsp',5);

    end
    methods
        function self = CTlungClass(varargin)
            self.initializeDat;
            if nargin > 2
                self.setcase(varargin{:});
            end
        end
        function setTag(self,varargin)
            N = numel(varargin)/2;
            if iscellstr(varargin) || isstring(varargin)
                if N == round(N)
                    for i = 1:N
                        self.dat.(varargin{2*i-1}).tag = varargin{2*i};
                    end
                else
                    warning('Invalid number of inputs.');
                end
            else
                warning('Invalid inputs.');
            end
        end
    end
end