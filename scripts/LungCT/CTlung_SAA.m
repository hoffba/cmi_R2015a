function res = CTlung_SAA(varargin)
% Optional inputs:
%   Nbins : (default = 10) number of slabs to average
%   endpad : (default = 15) % of volume on each end to ignore
%   dim : (default = 3) dimension of the map to plot
%   statstr : (default =  {'mean','std'}) statistics to calculate

dimstr = {'AP','LR','IS'};

[M,mask,Nbins,endpad,dim,statstr] = parseInputs(varargin);
N = size(M);
nstat = numel(statstr);

% Initialize results table
vnames = reshape(append(statstr,'_',{'grad','meanSlope'}'),1,[]);
vnames = reshape(append(dimstr(dim),'_',vnames'),1,[]);
res = table('Size',[1,numel(vnames)],...
    'VariableTypes',repmat({'double'},1,numel(vnames)),...
    'VariableNames',vnames);

for j = 1:numel(dim)

    % Calculate slabs for analysis
    slcdims = 1:3;
    slcdims(dim(j)) = [];
    ind = any(mask,slcdims);
    imin = find(ind,1);
    imax = find(ind,1,'last');
    d = imax-imin+1;
    imin = imin + floor(endpad*d/100);
    imax = imax - ceil(endpad*d/100);
    x = ceil(linspace(0,10,imax-imin+2));
    ind = zeros(N(dim(j)),1);
    ind(imin:imax) = x(2:end);

    % Generate slab statistics
    x = 1:Nbins;
    mvals = nan(Nbins,nstat);
    for i = x
        ii = ind==i;
        vals = getSliceVals(M,mask,dim(j),ii);
        for k = 1:nstat
            mvals(i,k) = feval(statstr{k},vals);
        end
    end

    % Polynomial fit over slab values
    for i = 1:nstat
        prestr = [dimstr{dim(j)},'_',statstr{i},'_'];
        % Gradient
        p = polyfit(x,mvals(:,i),1);
        res.([prestr,'grad']) = p(1);
        % Mean Slope
        res.([prestr,'meanSlope']) = (mvals(end,i) - mvals(1,i))/(Nbins-1);
    end

end


function vals = getSliceVals(M,mask,dim,ind)
switch dim
    case 1
        M = M(ind,:,:);
        mask = mask(ind,:,:);
    case 2
        M = M(:,ind,:);
        mask = mask(:,ind,:);
    case 3
        M = M(:,:,ind);
        mask = mask(:,:,ind);
end
vals = M(mask);

function [M,mask,Nbins,endpad,dim,statstr] = parseInputs(vargin)
p = inputParser;
addRequired(p,'M',@isnumeric);
addRequired(p,'mask',@islogical);
addOptional(p,'Nbins',10,@isnumeric);
addOptional(p,'endpad',15,@isnumeric);
addOptional(p,'dim',3,@isnumeric);
addOptional(p,'statstr',{'mean','std'},@iscellstr)
parse(p,vargin{:});
M = p.Results.M;
mask = p.Results.mask;
Nbins = p.Results.Nbins;
endpad = p.Results.endpad;
dim = p.Results.dim;
statstr = p.Results.statstr;


