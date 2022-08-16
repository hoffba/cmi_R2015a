classdef CTlungClass < handle
    properties (SetObservable, SetAccess=private, GetAccess=public)
        
        % Basic case info
        gui_flag = false;   % Switch for using GUI prompts
        swap_flag = false;  % Shows if images were automatically swapped at DICOM load
        id = "";        % Subject ID
        procdir = "";   % Processing directory
        elxdir = "";    % Elastix directory in procdir
        fn_log = "";    % File name of processing log
        
        opts = struct('quickreg',false,...  % Skip last resolution of warp
                      'swap_check',true,... % check for EXP/INS based on volume
                      'seg_method',4,...    % Segmentation method - 4=YACTA
                      'prm',[],...          % PRM settings
                      'tprm',struct('thresh',1:4,...
                                    'tmode','==',...
                                    'n',10*ones(1,3),...
                                    'gridsp',5));
        dat % Data structure containing relevant pipeline information:
            % i.e. filename tags, loaded image data, info, options
            % See method: initializeDat)
                
        ext = '.nii.gz';

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