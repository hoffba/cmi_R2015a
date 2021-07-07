classdef Mat3Dclass < handle
    % Class for general 3D matrix methods
    properties (SetAccess=protected, GetAccess=public, SetObservable)
        mat                 % Matrix of image/mask values
        dims = zeros(1,4);  % Dimensions of mat, 4-element vector
        voxsz = ones(1,3);  % Voxel dimensions
        voxsp = ones(1,3);  % Voxel spacing
        orient = [];        % Matrix spatial orientation affine matrix
        check = false;      % Check that matrix is available
    end
    
    methods
        function self = Mat3Dclass(obj)
            if nargin && isa(obj,'ImageClass')
                addlistener(obj,'voxsz','PostSet',@self.matchProps);
            end
        end
        function matchProps(self,src,evnt)
            if isprop(self,src.Name)
                self.(src.Name) = evnt.AffectedObject.(src.Name);
            end
        end
    end
    
    methods % Overloaded operators for Mat3Dclass and its subclasses
        % Overloaded +
        function self = plus(a,b)
            if isa(a,'Mat3Dclass')
                self = a;
            else
                self = b;
                b = a;
            end
            if isa(b,'MaskClass') && isa(self,'MaskClass')
                self.mat = self.mat & b.mat;
            elseif isa(b,'Mat3Dclass')
                self.mat = self.mat + b.mat;
            else
                self.mat = self.mat + b;
            end
        end
        % Overloaded -
        function self = minus(a,b)
            if isa(a,'Mat3Dclass')
                self = a;
            else
                self = -b;
                b = -a;
            end
            if isa(b,'MaskClass') && isa(self,'MaskClass')
                self.mat = self.mat & ~b.mat;
            elseif isa(b,'Mat3Dclass')
                self.mat = self.mat - b.mat;
            else
                self.mat = self.mat - b;
            end
        end
        % Overloaded -
        function self = uminus(self)
            if islogical(self.mat)
                self.mat = ~self.mat;
            else
                self.mat = -self.mat;
            end
        end
        % Overloaded .*
        function self = times(a,b)
            if isa(a,'Mat3Dclass')
                self = a;
            else
                self = b;
                b = a;
            end
            if isa(b,'MaskClass') && isa(self,'MaskClass')
                self.mat = self.mat & b.mat;
            elseif isa(b,'Mat3Dclass')
                self.mat = self.mat .* b.mat;
            else
                self.mat = self.mat .* b;
            end
        end
        % Overloaded ./
        function self = rdivide(a,b)
            if isa(b,'Mat3Dclass')
                self = b;
                b = b.mat;
            end
            if isa(a,'Mat3Dclass')
                self = a;
                a = a.mat;
            end
            if (islogical(a) && islogical(b))
                self.mat = a | b;
            else
                self.mat = a ./ b;
            end
        end
        % Overloaded /
        function self = mrdivide(a,b)
            % only for ease of use
            % does not perform matrix division on the image
            if isa(b,'Mat3Dclass')
                self = b;
                b = b.mat;
            end
            if isa(a,'Mat3Dclass')
                self = a;
                a = a.mat;
            end
            if (islogical(a) && islogical(b))
                self.mat = a | b;
            else
                self.mat = a ./ b;
            end
        end
        % Overloaded *
        function self = mtimes(a,b)
            % only for ease of use
            % does not perform matrix multiplication on the image
            if isa(a,'Mat3Dclass')
                self = a;
            else
                self = b;
                b = a;
            end
            if isa(b,'MaskClass') && isa(self,'MaskClass')
                self.mat = self.mat & b.mat;
            elseif isa(b,'Mat3Dclass')
                self.mat = self.mat .* b.mat;
            else
                self.mat = self.mat .* b;
            end
        end
        % Overloaded .^
        function self = power(a,b)
            if isa(b,'Mat3Dclass')
                self = b;
                b = b.mat;
            end
            if isa(a,'Mat3Dclass')
                self = a;
                a = a.mat;
            end
            if (islogical(a) && islogical(b))
                self.mat = a | b;
            else
                self.mat = a .^ b;
            end
        end
        % Overloaded ^
        function self = mpower(a,b)
            if isa(b,'Mat3Dclass')
                self = b;
                b = b.mat;
            end
            if isa(a,'Mat3Dclass')
                self = a;
                a = a.mat;
            end
            if (islogical(a) && islogical(b))
                self.mat = a | b;
            else
                self.mat = a .^ b;
            end
        end
        % Overloaded <
        function self = lt(a,b)
            if isa(b,'Mat3Dclass')
                self = b;
                b = b.mat;
            end
            if isa(a,'Mat3Dclass')
                self = a;
                a = a.mat;
            end
            self.mat = a < b;
        end
        % Overloaded >
        function self = gt(a,b)
            if isa(b,'Mat3Dclass')
                self = b;
                b = b.mat;
            end
            if isa(a,'Mat3Dclass')
                self = a;
                a = a.mat;
            end
            self.mat = a > b;
        end
        % Overloaded <=
        function self = le(a,b)
            if isa(b,'Mat3Dclass')
                self = b;
                b = b.mat;
            end
            if isa(a,'Mat3Dclass')
                self = a;
                a = a.mat;
            end
            self.mat = a <= b;
        end
        % Overloaded >=
        function self = ge(a,b)
            if isa(b,'Mat3Dclass')
                self = b;
                b = b.mat;
            end
            if isa(a,'Mat3Dclass')
                self = a;
                a = a.mat;
            end
            self.mat = a >= b;
        end
        % Overloaded ~=
        function self = ne(a,b)
            if isa(b,'Mat3Dclass')
                self = b;
                b = b.mat;
            end
            if isa(a,'Mat3Dclass')
                self = a;
                a = a.mat;
            end
            self.mat = a ~= b;
        end
        % Overloaded ==
        function self = eq(a,b)
            if isa(b,'Mat3Dclass')
                self = b;
                b = b.mat;
            end
            if isa(a,'Mat3Dclass')
                self = a;
                a = a.mat;
            end
            self.mat = a == b;
        end
        % Overloaded &
        function self = and(a,b)
            if isa(b,'Mat3Dclass')
                self = b;
                b = b.mat;
            end
            if isa(a,'Mat3Dclass')
                self = a;
                a = a.mat;
            end
            self.mat = a & b;
        end
        % Overloaded |
        function self = or(a,b)
            if isa(b,'Mat3Dclass')
                self = b;
                b = b.mat;
            end
            if isa(a,'Mat3Dclass')
                self = a;
                a = a.mat;
            end
            self.mat = a | b;
        end
        % Overloaded ~
        function self = not(self)
            self.mat = ~self.mat;
        end
        % Overloaded natural log
        function self = log(self)
            self.mat = log(self.mat);
        end
    end
end