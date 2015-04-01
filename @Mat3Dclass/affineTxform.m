% Mat3Dclass function
function Iout = affineTxform(self,varargin)
% Syntax: 
%   affineTxform
%       GUI input for transform options
%   affineTxform(A)
%       A : [4x4] transformation matrix
%   affineTxform('Name',Value,...)
%        ~~~'Name'~~~           ~~~Value~~~
%        'Rotate'               [1x3] vector of angles (deg)
%        'Translate'            [1x3] vector of displacements
%        'Scale'                [1x3] vector of scale factors
%        'Shear'                [1x6] vector of shear factors
%        'OutputLocations'      imref3d object with desired coordinates for output

% Validate inputs
Iout = [];
Cout = [];
A = [];
go = true;
if (nargin==1)
    A = varargin{1};
else % Need to create affine matrix
    if nargin==0
        answer = inputdlg({'Rotations (deg; X Y Z):','Translate (X Y Z):',...
                           'Scale (X Y Z):','Shear (YX ZX XY ZY XZ YZ):'},...
                          'Affine Transform:',1,...
                          {'0 0 0','0 0 0','1 1 1','0 0 0 0 0 0'});
        if ~isempty(answer)
            R = str2double(answer{1});
            T = str2double(answer{2});
            S = str2double(answer{3});
            K = str2double(answer{4});
        else
            go = false;
        end    
    elseif (mod(nargin,2)==0) && iscellstr(varargin(1:2:end))
        np = nargin/2;
        R = zeros(1,3); T = R;
        S = ones(1,3);
        K = zeros(1,6);
        for i = 1:np
            ii = 2*(i-1) + 1;
            switch varargin{ii}
                case 'Rotate'
                    R = varargin{ii+1};
                case 'Translate'
                    T = varargin{ii+1};
                case 'Scale'
                    S = varargin{ii+1};
                case 'Shear'
                    K = varargin{ii+1};
                case 'OutputLocations'
                    Cout = varargin{ii+1};
            end
        end
    else
        go = false;
    end
    if (length(R)==3) && (length(T)==3) && (length(S)==3) && (length(K)==6)
        A = genAffineTxform(T,R,S,K);
    end
end

% Perform the transformation
if go
    if ~isa(A,'affine3d')
        A = affine3d(A);
    end
    % Determine spatial coordinates for interpolation
    ext = ((self.dims(1:3)-1)/2 .* self.voxsz);
    ext = [-ext(:),ext(:)];
    Cin = imref3d(self.dims(1:3),ext(2,:),ext(1,:),ext(3,:));
    if ~isa(Cout,'imref3d')
        % If not input, assume same coordinates as object
        Cout = [];
    end
    % Transform the image:
    for i = 1:self.dims(4)
        if i==1
        end
        if isempty(Cout)
            % Allow image to grow
            Iout = imwarp(self.mat,Cin,A);
        else
            % Constrain image to defined size
            Iout = imwarp(self.mat,Cin,A,'OutputView',Cout);
        end
    end
end
