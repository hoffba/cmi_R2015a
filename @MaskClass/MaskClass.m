classdef MaskClass < Mat3Dclass
    properties (SetAccess=private, GetAccess=public)
        % Inherited properties
        %   mat         % Binary mask matrix
        %   dims        % Mask matrix dimensions
        %   check       % Check that matrix is available
    end
    methods
        function self = MaskClass(obj)
            if ~nargin
                obj = [];
            end
            self = self@Mat3Dclass(obj);
        end
    end
end