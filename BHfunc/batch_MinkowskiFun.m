%
% function batch_MinkowskiFun(fnout,fov,img,r,thresh,ind,varargin)
function batch_MinkowskiFun(fnout,img,mask,voxsz,ind,r,thresh,tmode,defVal,varargin)
% function batch_MinkowskiFun(fname,img,n,thresh,mask)
%   Runs a batch job to calculate local Minkowski Functionals over image
%   using defined image thresholds.
% Inputs:
%   fname
%   img
%   r
%   thresh
%   ind
%   tmode
%   defVal

if isempty(varargin)
    
    % This code starts the batch
	[~,fname] = fileparts(fnout);
    disp(['Starting batch Minkowski Functional analysis ... ',fname])
    batch(@batch_MinkowskiFun,0,{fnout,img,mask,voxsz,ind,r,thresh,tmode,defVal,1});
    disp(' ... done')
    
else
    
    % This code is run inside the batch
    disp(fnout);
    
    [MF,labels] = minkowskiFun(img,thresh,tmode,'mask',mask,'voxsz',voxsz,...
            'ind',ind,'n',r,'defVal',defVal);
    d = size(img);
    save(fnout,'MF','labels','thresh','ind','d','r','tmode','defVal','voxsz');
    
end