% ElxClass function
% Generate system call to elastix
function str = elxCmd(self,fname,mname,odir,varargin)
% Syntax: str = elxCmd(fname,mname,odir,'Name',Value,...)
%   * Automatically saves InitialTransform.txt and ElastixParameters.txt
%       determined by ElxClass properties.
% Inputs: fname = fixed image file name (.mhd)
%         mname = moving image file name (.mhd)
%         odir  = output directory
%         Name/Value pairs (optional)
%           'fMask' / FixedImageMask (.mhd)
%           'mMask' / MovingImageMask (.mhd)
%           'threads' / Number of threads used for optimization

str = '';
ns = length(self.Schedule);
if ns~=0
    mxCores = feature('numCores');
    % Parse variable inputs:
    p = inputParser;
    addRequired(p,'f',@ischar);
    addRequired(p,'m',@ischar);
    addRequired(p,'out',@ischar);
    addParameter(p,'fMask','',@ischar);
    addParameter(p,'mMask','',@ischar);
    addParameter(p,'ipp','',@ischar);
    addParameter(p,'threads',[],@(x)isnumeric(x)&&(x>0)&&(x<=mxCores));
    parse(p,fname,mname,odir,varargin{:})
    pp = p.Results;
    
    str = [fullfile(self.elxdir,'elastix'),...
           ' -f "',pp.f,'"',...
           ' -m "',pp.m,'"',...
           ' -out "',pp.out,'"'];
    
    % Add masks to string:
    mchk = false;
    if ~isempty(pp.fMask)
        str = [str,' -fMask "',pp.fMask,'"'];
        mchk = true;
    end
    if ~isempty(pp.mMask)
        str = [str,' -mMask "',pp.mMask,'"'];
    end
    
    if ~isempty(pp.ipp)
        str = [str,' -ipp "',pp.ipp,'"'];
    end
    
    % Optional threads input:
    if ~isempty(pp.threads)
        str = [str,' -threads ',num2str(pp.threads)];
    end
    
    % Save initial transform:
    if self.T0check && ~isempty(self.Tx0)
        pstr = self.saveTx0(fullfile(pp.out,'InitialTransform.txt'));
        str = [str,' -t0 "',pstr,'"'];
    end
    
    % Fix parameter structures
    for i = 1:ns
        if strcmp(self.getPar(i,'ImageSampler'),'RandomSparseMask') && ~mchk
            % Can't use RandomSparseMask if no mask was loaded.
            self.setPar(i,'ImageSampler','RandomCoordinate');
        end
    end
    
    % Save Parameter File(s):
    pstr = self.savePar(1:ns,fullfile(pp.out,'ElastixParameters.txt'));
    if ischar(pstr)
        pstr = {pstr};
    end
    str = [str,sprintf(' -p "%s"',pstr{:}),';'];
    
end

