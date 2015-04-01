% Run Elastix with given image files and list of optimization parameter files
function applyElastix(fixedFname,movingFname,parFname,varargin)

if nargin>2
    % First parse inputs:
    p = inputParser;
    addRequired(p,'FixedImage',@(x) (ischar(x) && (exist(x,'file')==2)));
    addRequired(p,'MovingImage',@(x) (ischar(x) && (exist(x,'file')==2)));
    addRequired(p,'Parameters',@iscellstr);
    addParamValue(p,'FixedMask','',@(x) (ischar(x) && (exist(x,'file')==2)));
    addParamValue(p,'MovingMask','',@(x) (ischar(x) && (exist(x,'file')==2)));
    addParamValue(p,'OutputDir','',@ischar);
    parse(p,fixedFname,movingFname,parFname,varargin{:});

    % Validate file types:
    [fpath,~,fext] = fileparts(p.Results.FixedImage);
    [~,~,mext] = fileparts(p.Results.MovingImage);
    nPar = length(parFname);
    pext = cell(1,nPar);
    for i = 1:length(parFname)
        [~,~,pext{i}] = fileparts(parFname{i});
    end
    if all(strcmp('.mhd',{fext,mext})) && all(strcmp('.txt',pext))
        % Generate call to Elastix:
        sysstr = ['/opt/elastix/bin/elastix',...
                  ' -f ',fixedFname,...
                  ' -m ',movingFname,...
                  sprintf(' -p %s',parFname{:})];
        if isempty(p.Results.OutputDir)
            sysstr = [sysstr,' -out ',fullfile(fpath,'elxreg')];
        else
            sysstr = [sysstr,' -out ',p.Results.OutputDir];
        end
        if ~isempty(p.Results.FixedMask)
            sysstr = [sysstr,' -fMask ',p.Results.FixedMask];
        end
        if ~isempty(p.Results.MovingMask)
            sysstr = [sysstr,' -mMask ',p.Results.MovingMask];
        end
        system(sysstr);
    end
end

