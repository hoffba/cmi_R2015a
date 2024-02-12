% Function for combining multiple transforms in Transformix
% Syntax: multiTform(tp)
% Inputs: tp = structure containing:
%               .fname  = string file name of TransformParameters.*.txt
%               .chain  = chain number (see below)
%               .im     = {nx4} cell array of image file names and output names
%                           * { NN , fname , outname , defVal }
%               .jac    = logical, to spit out SpatialJacobian
%         odir = directory for outputs
%         fmt = output image format
% 
function app = multiTform(tp,odir,fmt)

if nargin==0
    %% Load GUI
    app = multiTformGUI;
elseif isstruct(tp) && ~isempty(tp) && all(isfield(tp,{'fname','chain','im','jac'})) && ...
        ischar(odir) && exist(odir,'dir') && ischar(fmt) && ismember(fmt,{'.nii.gz','.nii','.mhd'})
    %% Run transforms
    t = tic;
    
    %% Determine which TP steps need [ nearest neighbor , linear ]:
    intrp = cell2mat(cellfun(@(x)[ (~isempty(x)&&any([x{:,1}])),...
                                   (isempty(x)||any(~[x{:,1}])) ] ,...
                             {tp(:).im}','UniformOutput',false));
    nx = find(intrp(:,1),1,'last');
    if ~isempty(nx)
        intrp(1:nx,1) = true;
    end
    nx = find(intrp(:,2),1,'last');
    if ~isempty(nx)
        intrp(1:nx,2) = true;
    end
    
    %% Initialize parameters
    logperm = 'wt';
    ncores = feature('NumCores');
    elxobj = ElxClass;
    elxobj.setFormat(fmt);
    logfname = fullfile(odir,'MultiTform.log');
    
    %% Loop over transform parameter steps, transform all images at each step
    ntp = length(tp);
    for i = 1:ntp
        %% Prep for new chain:
        chain = tp(i).chain;
        if (i==1) || (tp(i-1).chain~=chain)
            C = {};
            delete(fullfile(odir,'MultiTransPar.*.txt'));
            geochk = true;
        end
        
        %% Concatenate transform chain:
        nextname = tp(i).fname;
        bdir = fileparts(nextname);
        go = true;
        tC = {nextname};
        while go
            % Read TransformParameters file:
            disp(nextname);
            par = readTransformParam(nextname);
            nextname = par.InitialTransformParametersFileName;
            if strcmp(nextname,'NoInitialTransform')
                go = false;
            else
                [nextdir,tname,ext] = fileparts(nextname);
                % Make sure path is correct (copied folders/files)
                if ~exist(nextdir,'dir')
                    tname = fullfile(bdir,[tname,ext]);
                    if exist(tname,'file')
                        nextname = tname;
                    else
                        error('Could not find next transform file: %s',tname);
                    end
                end
                tC{end+1} = nextname;
            end
        end
        C = [C,fliplr(tC)];
        
        %% Read transformed image geometry from final transform
        if geochk % Do this once for each chain
            par = readTransformParam(C{1});
            sz = par.Size;
            sp = par.Spacing;
            affdir = par.Direction;
            orig = par.Origin;
            geochk = false;
        end
        
        %% Write chain to log file:
        jac = tp(i).jac;
        ncalls = max(jac,size(tp(i).im,1));
        if ncalls
            fid = fopen(logfname,logperm);
            if fid>2
                fprintf(fid,'\nImages:\n');
                if ~isempty(tp(i).im)
                    im = fliplr(tp(i).im(:,2:3))';
                    fprintf(fid,'   %s <-- %s\n',im{:});
                end
                if tp(i).jac
                    fprintf(fid,'   Jacobian\n');
                end
                fprintf(fid,'Transform Chain:\n');
                fprintf(fid,'   %s\n',C{:});
                fclose(fid);
            else
                error(['Could not open output directory: ',odr]);
            end
            logperm = 'at'; % change to append after first
        end
        
        %% Generate copies of all transforms in chain:
        ntot = length(C);
        for j = 1:ntot
            if j==ntot
                % Only need to set geometry in transformix input TP file:
                copytp(C{j},odir,j,intrp(i,:),sz,sp,affdir,orig,fmt);
            else
                copytp(C{j},odir,j,intrp(i,:));
            end
        end
        pause(0.1);
        
        %% Call transformix for selected images/jacobian:
        jac = tp(i).jac;
        for j = 1:ncalls
            chlink = sum([tp(1:i).chain]==chain);
            
            if ~isempty(tp(i).im)
                oname = tp(i).im{j,2};
                if strcmp(extractAfter(oname,'.'),'mhd')
                    % Fix MHD in case copied:
                    fixMHDnames(tp(i).im{j,2})
                end
                if (size(tp(i).im,2)==1) || isempty(tp(i).im{j,3})
                    outfn = fullfile(odir,sprintf('TransformedImage.%02u.%02u.%02u.mhd',chain,chlink,j));
                else
                    outfn = fullfile(odir,tp(i).im{j,3});
                end
                inputs = {'in',tp(i).im{j,2},'outfn',outfn};
                nn = tp(i).im{j,1};
                [~,str,~] = fileparts(outfn);
            else
                nn = false;
                str = 'spatialJacobian';
                inputs = {};
            end
            
            tpcall = tpfname(odir,ntot,nn);
            setDefVal(tpcall,tp(i).im{j,4});
            
            fprintf('   Transforming: %s\n',str);
            
            str = elxobj.sysCmd(odir,inputs{:},'tp',tpcall,'jac',jac,...
                'wait',true,'threads',ncores);
            system(str);
            
            if jac
                jac = false; % only calculate once per transform
                fnames = dir(fullfile(odir,'spatialJacobian.*'));
                for ifn = 1:length(fnames)
                    fnout = insertBefore(fnames(i).name,'.',sprintf('.%02u.%02u',chain,chlink));
                    movefile(fullfile(odir,fnames(ifn).name),fullfile(odir,fnout));
                end
            end
            
            % Clear temporary files:
            fnames = dir(fullfile(odir,'result.*'));
            if ~isempty(fnames)
                fnames = fullfile(odir,{fnames.name});
                delete(fnames{:});
            end

            % Copy log file:
            movefile(fullfile(odir,'transformix.log'),...
                fullfile(odir,sprintf('transformix.%02u.%02u.%02u.log',...
                chain,chlink,j)));
        end
        
    end
    fixMHDnames(odir);
    h = toc(t);
end

function str = tpfname(odir,i,nn)
    nnstr = '';
    if nn
        nnstr = '_NN';
    end
    str = fullfile(odir,sprintf('MultiTransPar.%02u%s.txt',i,nnstr));

function copytp(fname,odir,i,nn,sz,sp,affdir,orig,fmt)
% Copies transform parameter file to new directory and adjust relevant parameters:
% Inputs:
%   fname   = transform parameter file to copy
%   odir    = directory to copy file to
%   i       = step in transform chain
%   fin     = T/F indicating last transform in chain
%   nn      = T/F nearest neighbor interpolation
%   sz      = interpolating matrix size
%   sp      = interpolating image voxel spacing
%   orig    = interpolating image geometric origin
%   fmt     = image output format (i.e. 'nii' or 'mhd')

    % Read file:
    p = readTransformParam(fname);

    if isfield(p,'HowToCombineTransforms') && strcmp(p.HowToCombineTransforms(1),'c')
        p.HowToCombineTransforms = 'Compose'; % Must be capitalized
    end

    % Parameters for final transform in chain
    if nargin==9
        p.Size = sz;
        p.Spacing = sp;
        p.Direction = affdir;
        p.Origin = orig;
        if endsWith(fmt,'.gz')
            p.CompressResultImage = 'true';
            fmt(end-2:end) = [];
        end
        p.ResultImageFormat = fmt(2:end);
    end

    if i==1
        p.InitialTransformParametersFileName = 'NoInitialTransform';
    end 
    
    if nn(1) % Nearest Neighbor flag
        if i>1
            p.InitialTransformParametersFileName = tpfname(odir,i-1,true);
        end
        p.ResampleInterpolator = 'FinalNearestNeighborInterpolator';
        writeTransformParam(tpfname(odir,i,true),p);
    end
    
    if nn(2) % Linear flag
        if i>1
            p.InitialTransformParametersFileName = tpfname(odir,i-1,false);
        end
        p.ResampleInterpolator = 'FinalBSplineInterpolator';
        writeTransformParam(tpfname(odir,i,false),p);
    end
    
% Set DefaultPixelValue for existing parameter file:
function setDefVal(fname,defVal)
    % Read current parameter file:
    fid = fopen(fname,'rt');str = fread(fid,'*char')';fclose(fid);
    % Change DefaultPixelValue to desired value:
    str = regexprep(str,'\(DefaultPixelValue [^)]*',...
        ['\(DefaultPixelValue ',sprintf('%.6f',defVal)]);
    % Re-write file:
    fid = fopen(fname,'wt');fwrite(fid,str,'char');fclose(fid);
