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
    addRequired(p,'odir',@ischar);
    addParameter(p,'fMask','',@ischar);
    addParameter(p,'mMask','',@ischar);
    addParameter(p,'ipp','',@ischar);
    addParameter(p,'threads',[],@(x)isnumeric(x)&&(x>0)&&(x<=mxCores));
    parse(p,fname,mname,odir,varargin{:})
    pp = p.Results;
    
    str = [fullfile(self.elxdir,'elastix'),...
           ' -f "',pp.f,'"',...
           ' -m "',pp.m,'"',...
           ' -out "',pp.odir,'"'];
    
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
        pstr = self.saveTx0(fullfile(pp.odir,'InitialTransform.txt'));
        str = [str,' -t0 "',pstr,'"'];
    end
    
    
    
    
    
    
    
    % Copy and adjust secondary initial transform guess:
    if ischar(self.Tx0guess) && exist(self.Tx0guess,'file')
        fname = self.Tx0guess;
        it = true;
        ct = 0;
        while it
            
            % Read file:
            fid = fopen(fname,'rt');
            str = fread(fid,ind,'*char')';
            fclose(fid);
            
            % Find sub-initial transform:
            [ext,tok] = regexp(str,'\(InitialTransformParametersFileName \"(.*?)\"\)',...
                'tokenExtents','tokens');
            fname = tok{1}{1};
            it = ~strcmp(fname,'NoInitialTransform');
            if it 
                if exist(fname,'file')
                    str = [ str(1:ext{1}(1)-1) ,...
                            fullfile(pp.odir,sprintf('InitialGuess_%u.txt',ct+1)) ,...
                            str(ext{1}(2)+1:end) ];
                else
                    error(['Missing transform file: ',fname]);
                end
            end
            
            % Save new file to output directory:
            fid = fopen(fullfile(self.odir,sprintf('InitialGuess_%u.txt',ct)),'wt');
            fwrite(fid,str,'char');
            fclose(fid);
        end
    end
    
    
    
    
    
    
    
    
    
    % Fix parameter structures
    for i = 1:ns
        if strcmp(self.getPar(i,'ImageSampler'),'RandomSparseMask') && ~mchk
            % Can't use RandomSparseMask if no mask was loaded.
            self.setPar(i,'ImageSampler','RandomCoordinate');
        end
    end
    
    % Save Parameter File(s):
    pstr = self.savePar(1:ns,fullfile(pp.odir,'ElastixParameters.txt'));
    if ischar(pstr)
        pstr = {pstr};
    end
    str = [str,sprintf(' -p "%s"',pstr{:})];
    
end

