% ImageClass function
function [MFout,labels] = calcMF(self,vec,thresh,n,varargin)
% Calculate Minkowski Functionals
% Inputs:   vec = 4D image index to perform analysis on
%           optional Name/Value pairs:
%               'Thresh'  = vector of thresholds to use
%               'Window'  = radius (voxels) of moving window
%                               *(use inf if you want whole-image)
%               'ApplyMask' = check for use of existing VOI
%               'ImgCheck'  = check to analyze image (true) or VOI (false)

MFout = [];
labels = {};

p = inputParser;
addRequired(p,'Vec',@isscalar);
addRequired(p,'Thresh',@isnumeric);
addRequired(p,'Window',@isvector);
addParameter(p,'ApplyMask',true,@islogical);
parse(p,vec,thresh,n,varargin{:});
pp = p.Results;

varargin = {};
stat = false;
if (pp.Vec==0)
    if self.mask.check
        % Use existing VOI
        timg = double(self.mask.mat);
        pp.Thresh = 0.5;
        stat = true;
    else
        warning('No VOI exists for Minkowski Functional analysis.');
    end
elseif self.check
    timg = self.mat(:,:,:,pp.Vec);
    stat = true;
    if pp.ApplyMask
        varargin = {self.mask.mat};
    end
end

if stat
    disp('Starting batch MF analysis ...');
    batch(@minkowskiFun,2,[{timg,pp.Window,pp.Thresh},varargin]);
    disp(' ... done')
%     [MFout,labels] = minkowskiFun(timg,pp.Window,pp.Thresh,varargin{:});
%     if ~isempty(MFout) && (length(size(MFout))==5)
%         d = size(MFout);
%         labels = strcat(repmat(labels,d(4),1),...
%                         repmat(cellfun(@num2str,num2cell(1:d(4))',...
%                                        'UniformOutput',false),1,d(5)));
%         self.imgAppend(reshape(MFout,[d(1:3),prod(d(4:end))]),labels);
%     end
end
