%
% function batch_MinkowskiFun(fnout,fov,img,r,thresh,ind,varargin)
function batch_MinkowskiFun(fnout,img,varargin)
% function batch_MinkowskiFun(fname,img,n,thresh,mask)
%   Runs a batch job to calculate local Minkowski Functionals over image
%   using defined image thresholds.
% Inputs:
%   fname
%   img
%   Name/Value pairs for MF options (see minkowskiFun.m)

if (nargin/2)==round(nargin/2)
    
    % This code starts the batch
	[~,fname] = fileparts(fnout);
    disp(['Starting batch Minkowski Functional analysis ... ',fname])
    batch(@batch_MinkowskiFun,0,[{fnout,img},varargin,{1}]);
    disp(' ... done')
    
else
    
    % This code is run inside the batch
    disp(fnout);
    
    [MF,p] = minkowskiFun(img,varargin{:});
    save(fnout,'MF','p');
    
end