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
    pstr = '';
    if self.T0check && ~isempty(self.Tx0)
        pstr = self.saveTx0(fullfile(pp.odir,'InitialTransform.txt'));
    end
    
    % Copy and adjust secondary initial transform guess:
    if ~isempty(self.Tx0guess)
        C = struct2cell(self.Tx0guess);
        i = [C{1,:}];
        fnames = cellfun(@(x,y)fullfile(x,y),C(3,i),C(2,i),'UniformOutput',false)';
        nf = length(fnames);
        for i = 1:nf
            
            % Read parameter file:
            if exist(fnames{i},'file')
                fid = fopen(fnames{i},'rt');
                txt = fread(fid,inf,'*char')';
                fclose(fid);
            else
                error(['Missing transform file: ',fname]);
            end
            
            % Change initial transform file name:
            ext = regexp(txt,'\(InitialTransformParametersFileName \"(.*?)\"\)','tokenExtents');
            if (i==nf)
                if isempty(itname)
                    % No initial transform
                    igname = 'NoInitialTransform';
                else
                    igname = itname;
                end
            else
                igname = fullfile(pp.odir,sprintf('InitialGuess_%u.txt',i+1));
            end
            txt = [ txt(1:ext{1}(1)-1) , igname , txt(ext{1}(2)+1:end) ];
            
            % Save copy into new output directory for this optimization:
            igname = fullfile(pp.odir,sprintf('InitialGuess_%u.txt',i));
            if i==1
                itname = pstr;
                pstr = igname;
            end
            fid = fopen(igname,'wt');
            fwrite(fid,txt,'char');
            fclose(fid);
        end
    end
    
    % Set initial transform for this elastix call:
    if ~isempty(pstr)
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
    pstr = self.savePar(1:ns,fullfile(pp.odir,'ElastixParameters.txt'));
    if ischar(pstr)
        pstr = {pstr};
    end
    str = [str,sprintf(' -p "%s"',pstr{:})];
    
end

