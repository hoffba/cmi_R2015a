function varargout = cmi_gradient(f,varargin)
%GRADIENT Approximate gradient.
%   [FX,FY] = GRADIENT(F) returns the numerical gradient of the
%   matrix F. FX corresponds to dF/dx, the differences in x (horizontal) 
%   direction. FY corresponds to dF/dy, the differences in y (vertical) 
%   direction. The spacing between points in each direction is assumed to 
%   be one. When F is a vector, DF = GRADIENT(F)is the 1-D gradient.
%
%   [FX,FY] = GRADIENT(F,H), where H is a scalar, uses H as the
%   spacing between points in each direction.
%
%   [FX,FY] = GRADIENT(F,HX,HY), when F is 2-D, uses the spacing
%   specified by HX and HY. HX and HY can either be scalars to specify
%   the spacing between coordinates or vectors to specify the
%   coordinates of the points.  If HX and HY are vectors, their length
%   must match the corresponding dimension of F.
%
%   [FX,FY,FZ] = GRADIENT(F), when F is a 3-D array, returns the
%   numerical gradient of F. FZ corresponds to dF/dz, the differences
%   in the z direction. GRADIENT(F,H), where H is a scalar, 
%   uses H as the spacing between points in each direction.
%
%   [FX,FY,FZ] = GRADIENT(F,HX,HY,HZ) uses the spacing given by
%   HX, HY, HZ. 
%
%   [FX,FY,FZ,...] = GRADIENT(F,...) extends similarly when F is N-D
%   and must be invoked with N outputs and either 2 or N+1 inputs.
%
%   Note: The first output FX is always the gradient along the 2nd
%   dimension of F, going across columns.  The second output FY is always
%   the gradient along the 1st dimension of F, going across rows.  For the
%   third output FZ and the outputs that follow, the Nth output is the
%   gradient along the Nth dimension of F.
%
%   Examples:
%       [x,y] = meshgrid(-2:.2:2, -2:.2:2);
%       z = x .* exp(-x.^2 - y.^2);
%       [px,py] = gradient(z,.2,.2);
%       contour(z), hold on, quiver(px,py), hold off
%
%   Class support for input F:
%      float: double, single
%
%   See also DIFF, DEL2.

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 5.17.4.8 $  $Date: 2011/05/17 02:22:19 $

[err,f,ndim,loc,rflag] = parse_inputs(f,varargin);
if err, error(message('MATLAB:gradient:InvalidInputs')); end

% Loop over each dimension. Permute so that the gradient is always taken along
% the columns.

if ndim == 1
  perm = [1 2];
else
  perm = [2:ndim 1]; % Cyclic permutation
end

for k = 1:ndim
   [n,p] = size(f);
   h = loc{k}(:); 
   vchk = length(h)>1;
   g  = zeros(size(f),class(f)); % case of singleton dimension

   % Take forward differences on left and right edges
   if n > 1
      if vchk
          d = [h(2)-h(1) , h(end)-h(end-1)];
      else
          d = h*ones(1,2);
      end
      g(1,:) = (f(2,:) - f(1,:))/d(1);
      g(n,:) = (f(n,:) - f(n-1,:))/d(2);
   end

   % Take centered differences on interior points
   if n > 2
      if vchk
          d = h(3:end) - h(1:end-2);
          g(2:n-1,:) = (f(3:end,:)-f(1:end-2,:))./d(:,ones(p,1));
      else
          d = 2*h;
          g(2,:) = (f(3,:)-f(1,:))/d;
          g(n-1,:) = (f(end,:)-f(end-2,:))/d;
      end
   end
   
   % Take five-point stencil approximation for interior points
   if ~vchk && (n > 4)
      d = 12*h;
      g(3:n-2,:) = (f(1:n-4,:)-8*f(2:n-3,:)+8*f(4:n-1,:)-f(5:n,:))/d;
   end

   varargout{k} = ipermute(g,[k:max(ndim,2) 1:k-1]);

   % Set up for next pass through the loop
   f = permute(f,perm);
end 

% Swap 1 and 2 since x is the second dimension and y is the first.
if ndim>1
  tmp = varargout{1};
  varargout{1} = varargout{2};
  varargout{2} = tmp;
end

if rflag, varargout{1} = varargout{1}.'; end


%-------------------------------------------------------
function [err,f,ndim,loc,rflag] = parse_inputs(f,v)
%PARSE_INPUTS
%   [ERR,F,LOC,RFLAG] = PARSE_INPUTS(F,V) returns the spacing
%   LOC along the x,y,z,... directions and a row vector
%   flag RFLAG. ERR will be true if there is an error.

err = false;
loc = {};
nin = length(v)+1;

% Flag vector case and row vector case.
ndim = ndims(f);
vflag = 0; rflag = 0;
if iscolumn(f)
   ndim = 1; vflag = 1; 
elseif isrow(f) % Treat row vector as a column vector
   ndim = 1; vflag = 1; rflag = 1;
   f = f.';
end;

% Default step sizes: hx = hy = hz = 1
if nin == 1, % gradient(f)
   loc = cell(1, ndims(f));
   for k = 1:ndims(f)
      loc(k) = {1};
   end;

elseif (nin == 2) % gradient(f,h)
   % Expand scalar step size
   if (length(v{1})==1)
      loc = cell(1, ndims(f)); 
      for k = 1:ndims(f)
         h = v{1};
         loc(k) = {h};
      end;
   % Check for vector case
   elseif vflag
      loc(1) = v(1);
   else
      err = true;
   end

elseif ndims(f) == numel(v), % gradient(f,hx,hy,hz,...)
   % Swap 1 and 2 since x is the second dimension and y is the first.
   loc = v;
   if ndim>1
     tmp = loc{1};
     loc{1} = loc{2};
     loc{2} = tmp;
   end

   % replace any scalar step-size with corresponding position vector
   for k = 1:ndims(f)
      if length(loc{k})==1
         loc{k} = loc{k};
      end;
   end;

else
   err = true;

end
