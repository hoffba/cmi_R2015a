function func = map3d_select(varargin)
% Function for selecting and running a moving window process
% Inputs follow Name/Value syntax:
%       Name = function name, i.e. 'kurtosis' OR 'MF'
%       Value = cell array of additional inputs for function handle
%           i.e. for 'MF': {voxsz,ord}

guichk = nargin==0;
opts = {'kurtosis', { @(img,mask)kurtosis(img(mask)) };...
        'perc',     { @(img,mask,x)prctile(img(mask),x), 15 };... % default to 15
        'mean',     { @(img,mask)mean(img(mask)) };...
        'mean-all', { @(img,mask)mean(img,'all') };...
        'std',      { @(img,mask)std(img(mask)) };...
        'std-all',  { @(img,mask)std(img,0,'all') };...
        'MF',       { @(img,mask,voxsz,ord)calcMF3D(img,voxsz,ord), ones(1,3), 0:3 }};

if guichk
    % Use GUI to select function handle
    answer = listdlg('ListString',opts(:,1));
    answer = opts(answer,1);
    [~,Locb] = ismember(answer,opts(:,1));
    func = opts(Locb,2);
elseif nargin/2 == round(nargin/2) % check for even number of inputs
    answer = varargin(1:2:end);
    inputs = varargin(2:2:end);
    [Lia,Locb] = ismember(answer,opts(:,1));
    nf = numel(answer);
    func = cell(nf,1);
    for i = 1:nf
        if ~Lia(i)
            warning('Could not find match for function: %s',answer{i}); 
        elseif ~isempty(inputs{i}) && (numel(inputs{i}) ~= (numel(opts{Locb(i),2})-1))
            warning('Number of inputs for function %s is invalid',answer{i});
        else
            func{i} = opts{Locb(i),2};
            if ~isempty(inputs{i})
                func{i}(2:end) = inputs{i};
            end
        end
    end
else
    warning('Invalid number of inputs.');
    return;
end
