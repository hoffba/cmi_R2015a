classdef ElxClass < handle
    properties (SetObservable, SetAccess=private, GetAccess=public)
        
        elxdir = '/opt/elastix/bin/';
        
        % Stack of multi-parameter optimizations
        T0check = false; % Determines whether to apply initial transform during Elastix call
        Tx0            % Initial transformation (before optimization)
        Schedule = {}; % Schedule of coregistration parameters
        
    end
    methods
        % Constructor
        function self = ElxClass
            if ismac
                % You won't need the full path if you've installed Elastix
                % correctly
                self.elxdir = '';
            end
        end
    end
end