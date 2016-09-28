% Function for combining multiple transforms in Transformix
% Syntax: multiTform(tp)
% Inputs: tp = structure containging:
%               .fname  = string file name of TransformParameters.*.txt
%               .chain  = chain number (see below)
%               .im     = {nx3} cell array of image file names and output names
%                           * { NN , fname , outname }
%               .jac    = logical, to spit out SpatialJacobian
%         odir = directory for outputs
% 
function h = multiTform(tp,odir)

if nargin==0
    h = multiTformGUI;
elseif isstruct(tp) && all(isfield(tp,{'fname','chain','im','jac'})) && ...
        ~isempty(tp) && ischar(odir) && exist(odir,'dir')
    t = tic;
    
    % Determine which TP steps need nearest neighbor / linear:
    intrp = cell2mat(cellfun(@(x)[ (~isempty(x)&&any([x{:,1}])),...
                                   (~isempty(x)&&any(~[x{:,1}])) ] ,...
                             {tp(:).im}','UniformOutput',false));
    nx = find(intrp(:,1),1,'last');
    if ~isempty(nx)
        intrp(1:nx,1) = true;
    end
    nx = find(intrp(:,2),1,'last');
    if ~isempty(nx)
        intrp(1:nx,2) = true;
    end
    
    logperm = 'wt';
    ncores = feature('NumCores');
    elxobj = ElxClass;
    logfname = fullfile(odir,'MultiTform.log');
    
    ntp = length(tp);
    C = {};
    % Loop over transform parameter steps, transform all images at each step
    for i = 1:ntp
        
        % Concatenate transform chain:
        nextname = tp(i).fname;
        bdir = fileparts(nextname);
        go = true;
        tC = {nextname};
        while go
            % Read txt file:
            disp(nextname);
            fid = fopen(nextname,'r');
            if fid>2
                str = fread(fid,'*char')';
                fclose(fid);
            else
                error(['Could not open file: ',nextname]);
            end
            % Search for initial transform:
            str = regexp(str,'\(InitialTransformParametersFileName "([^"]*)','tokens');
            nextname = str{1}{1};
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
        
        % Read transformed image geometry for interpolation:
        if i==1
            fid = fopen(C{1},'r');
            sz = []; sp = []; orig = [];
            while isempty(sz) || isempty(sp) || isempty(orig)
                str = fgetl(fid);
                tok = regexp(str,'\((Size|Spacing|Origin) ([^)]*)','tokens');
                if ~isempty(tok)
                    switch tok{1}{1}
                        case 'Size'
                            sz = str2num(tok{1}{2});
                        case 'Spacing'
                            sp = str2num(tok{1}{2});
                        case 'Origin'
                            orig = str2num(tok{1}{2});
                    end
                end
            end
            fclose(fid);
        end
        
        % Write to log file:
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
        
        % Generate copies of all transforms in chain:
        ntot = length(C);
        for j = 1:ntot
            if j==ntot
                % Only need to set in transformix input TP file:
                copytp(C{j},odir,j,intrp(i,:),sz,sp,orig);
            else
                copytp(C{j},odir,j,intrp(i,:));
            end
        end
        pause(0.1);
        
        % Call transformix for selected images/jacobian:
        jac = tp(i).jac;
        for j = 1:ncalls
            
            chain = tp(i).chain;
            chlink = sum([tp(1:i).chain]==tp(i).chain);
            
            if ~isempty(tp(i).im)
                % Fix MHD in case copied:
                fixMHDnames(tp(i).im{j,2})
                if (size(tp(i).im,2)==1) || isempty(tp(i).im{j,3})
                    fmt = 'TransformedImage.%02u.%02u.%02u.mhd';
                    outfn = fullfile(odir,sprintf(fmt,chain,chlink,j));
                else
                    outfn = fullfile(odir,tp(i).im{j,3});
                end
                inputs = {'in',tp(i).im{j,2},'outfn',outfn};
            end
            
            tpcall = tpfname(odir,ntot,tp(i).im{j,1});
            [~,str,~] = fileparts(outfn);
            fprintf('   Transforming: %s\n',str);
            
            str = elxobj.sysCmd(odir,inputs{:},'tp',tpcall,'jac',jac,...
                'wait',true,'threads',ncores);
            system(str);
            
            if jac
                jac = false; % only calculate once per transform
                movefile(fullfile(odir,'spatialJacobian.mhd'),...
                    fullfile(odir,sprintf('spatialJacobian.%02u.%02u.mhd',...
                    chain,chlink)));
                movefile(fullfile(odir,'spatialJacobian.raw'),...
                    fullfile(odir,sprintf('spatialJacobian.%02u.%02u.raw',...
                    chain,chlink)));
            end
            % Clear temporary files:
            if exist(fullfile(odir,'result.mhd'),'file')
                delete(fullfile(odir,'result.mhd'),fullfile(odir,'result.raw'));
            end
            % Copy log file:
            movefile(fullfile(odir,'transformix.log'),...
                fullfile(odir,sprintf('transformix.%02u.%02u.%02u.log',...
                chain,chlink,j)));
        end
        
        % Prep for new chain:
        if (i==ntp) || (tp(i+1).chain~=tp(i).chain)
            C = {};
            delete(fullfile(odir,'MultiTransPar.*.txt'));
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

function copytp(fname,odir,i,nn,sz,sp,orig)
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

% Read file:
fid = fopen(fname,'r'); str = fread(fid,'*char')'; fclose(fid);
% Temporary fix for previous coregistrations:
% * (found that HowToCombineTransforms is case sensitive)
str = regexprep(str,'\(HowToCombineTransforms ([^)]*)',...
    '\(HowToCombineTransforms "Compose"');
if nargin==7 % Only needs to be done on parameter file that is called
    str = regexprep(str,'\(Size ([^)]*)',...
        ['\(Size ',sprintf('%u %u %u',sz)]);
    str = regexprep(str,'\(Spacing ([^)]*)',...
        ['\(Spacing ',sprintf('%.6f %.6f %.6f',sp)]);
    str = regexprep(str,'\(Origin ([^)]*)',...
        ['\(Origin ',sprintf('%.10f %.10f %.10f',orig)]);
end
if nn(1) % Nearest Neighbor flag
    if i==1
        itname = 'NoInitialTransform';
    else
        itname = tpfname(odir,i-1,true);
    end
    str = regexprep(str,'\(InitialTransformParametersFileName "([^"]*)',...
        ['\(InitialTransformParametersFileName "',itname]);
    str = regexprep(str,'\(ResampleInterpolator ([^)]*)',...
        '\(ResampleInterpolator "FinalNearestNeighborInterpolator"');
    fid = fopen(tpfname(odir,i,true),'w');
    fwrite(fid,str,'char');
    fclose(fid);
end
if nn(2) % Linear flag
    if i==1
        itname = 'NoInitialTransform';
    else
        itname = tpfname(odir,i-1,false);
    end
    str = regexprep(str,'\(InitialTransformParametersFileName "([^"]*)',...
        ['\(InitialTransformParametersFileName "',itname]);
    str = regexprep(str,'\(ResampleInterpolator ([^)]*)',...
        '\(ResampleInterpolator "FinalBSplineInterpolator"');
    fid = fopen(tpfname(odir,i,false),'w');
    fwrite(fid,str,'char');
    fclose(fid);
end
    