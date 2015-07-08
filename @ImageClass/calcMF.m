% ImageClass function
function [MF,labels] = calcMF(self,vec,thresh,n,varargin)
% Calculate Minkowski Functionals
% Inputs:   vec = 4D image index to perform analysis on
%           ithresh = vector of thresholds to use
%           n = window radius (voxels) (inf = global analysis)
%           optional Name/Value pairs:
%               'ApplyMask' = check for use of existing VOI
%               'tmode'     = string describing thresholding operator
%                               - starts with 'p' --> percentiles
%                               - {'<','<=','>','>=','==','~='}
%               'defVal'    = logical value of voxels outside the VOI

MF = [];
labels = {};

p = inputParser;
addRequired(p,'vec',@(x)isscalar(x)&&(x>0));
addRequired(p,'thresh',@isnumeric);
addRequired(p,'n',@isvector);
addParameter(p,'ApplyMask',true,@islogical);
addParameter(p,'tmode','>',@ischar);
addParameter(p,'defVal',false,@islogical);

parse(p,vec,thresh,n,varargin{:});
pp = p.Results;

if self.check
    if pp.vec
        timg = self.mat(:,:,:,pp.vec);
    else % Perform analysis on VOI
        timg = self.mask.mat;
        pp.ApplyMask = false;
        pp.thresh = 0.5;
        pp.tmode = '>';
    end
    if pp.ApplyMask && self.mask.check
        mask = self.mask.mat;
    else
        mask = [];
    end
    gchk = any(isinf(pp.n));
    
    if gchk % Global analysis:

        % Calculate MF:
        [MF,labels] = minkowskiFun(timg,pp.thresh,pp.tmode,...
            'n',pp.n,'voxsz',self.voxsz,'mask',mask,'defVal',pp.defVal,'prog',true);
            
        if nargout==0
            % Save results to base workspace:
            pp.MF = MF;
            pp.labels = labels;
            assignin('base','MFvals',pp);
        end
        
    else
        
        % Run in batch - will take too long to wait for
        
        fnout = fullfile(self.dir,[self.name,'_MF',datestr(now(),'yyyymmddHHMM'),'.mat']);

        % Generate grid indices:
        mchk = isempty(mask);
        if mchk
            mmin = [ find(max(max(mask,[],2),[],3),1) ,...
                     find(max(max(mask,[],1),[],3),1) ,...
                     find(max(max(mask,[],1),[],2),1) ];
        else
            mmin = ones(1,3);
        end
        gridmask = false(self.dims(1:3));
        gridmask( mmin(1):pp.n(1):end , ...
                  mmin(2):pp.n(2):end , ...
                  mmin(3):pp.n(3):end ) = true;
        if mchk
            gridmask = gridmask & mask;
        end
        ind = find(gridmask);
        
        batch_MinkowskiFun(fnout,timg,mask,self.voxsz,ind,pp.n,pp.thresh,pp.tmode,pp.defVal);
    end
end

