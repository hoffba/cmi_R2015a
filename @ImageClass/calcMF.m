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
addRequired(p,'Vec',@(x)isscalar(x)&&(x>0));
addRequired(p,'Thresh',@isnumeric);
addRequired(p,'Window',@isvector);
addParameter(p,'ApplyMask',true,@islogical);
parse(p,vec,thresh,n,varargin{:});
pp = p.Results;
r = pp.Window;

if self.check
    
    timg = self.mat(:,:,:,pp.Vec);
    gchk = any(isinf(r));
    mchk = pp.ApplyMask && self.mask.check;
    
    if gchk
        
        nth = length(pp.Thresh);
        if length(r) == 2
            nmf = 3;
        else
            nmf = 4;
        end

        % Initialize output matrix:
        MFout = zeros([nth,nmf]);

        % Loop over desired image thresholds:
        hw = waitbar(0,'Calculating Gobal Minkowski Functionals ...');
        for ith = 1:nth
            BW = timg > pp.Thresh(ith);
            if mchk
                BW = BW & self.mask.mat;
            end
            [MFout(ith,:),labels] = minkowskiFun(BW,r,[],self.voxsz);
            waitbar(ith/nth,hw,['Calculating Gobal Minkowski Functionals ... ',...
                                num2str(pp.Thresh(ith)),' done']);
        end
        delete(hw);

        if nargout==0
            % Display global results:
            vals = [{'Thresh'},num2cell(pp.Thresh);...
                    labels(:),num2cell(squeeze(MFout'))];
            assignin('base','MRvals',vals);
            msgbox(sprintf(['%10s:',repmat(' %8f',1,nth)],vals{:}),'Minkowski Functionals');
        end
    else
        
        % Run in batch - will take too long to wait for
        
        % User decides where to save results:
        [fname,fpath] = uiputfile('*.mat','Save MF Results',self.name);
        
        ind = [];
        if mchk
            mask = self.mask.mat;
            mmin = [ find(max(max(mask,[],2),[],3),1) ,...
                     find(max(max(mask,[],1),[],3),1) ,...
                     find(max(max(mask,[],1),[],2),1) ];
            gridmask = false(self.dims(1:3));
            gridmask( mmin(1):r(1):end , ...
                      mmin(2):r(2):end , ...
                      mmin(3):r(3):end ) = true;
            ind = find(self.mask.mat & gridmask);
        else
            mask = true(self.mask.dims(1:3));
        end
        batch_MinkowskiFun(fullfile(fpath,fname),timg,mask,self.voxsz,ind,...
                            r,pp.Thresh);
        
    end
end

