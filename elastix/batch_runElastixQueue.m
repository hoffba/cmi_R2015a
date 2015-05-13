
function batch_runElastixQueue(fname,varargin)
% Start / run batch Elastix queue from specified file:
% Inputs:   fname = File name of queue

if nargin<5
    % This code starts the batch
    disp('Starting batch heterogeneity analysis process ...')
    batch(@batch_runElastixQueue,0,{fname,1});
    disp(' ... done')
elseif varargin{1}==1
    
    % Start next process:
    
    
end
