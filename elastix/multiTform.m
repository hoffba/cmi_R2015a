% Function for combining multiple transforms in Transformix
% Syntax: multiTform(tp)
% Inputs: tp = structure containging:
%               .fname  = string file name of TransformParameters.*.txt
%               .chain  = chain number (see below)
%               .im     = cell array of image file names
%               .jac    = logical, to spit out SpatialJacobian
%         odir = directory for outputs
% 
function h = multiTform(tp,odir)

if nargin==0
    h = multiTformGUI;
elseif isstruct(tp) && all(isfield(tp,{'fname','chain','im','jac'})) && ...
        ~isempty(tp) && ischar(odir) && exist(odir,'dir')
    t = tic;
    ncores = feature('NumCores');
    elxobj = ElxClass;
    logfname = fullfile(odir,'MultiTform.log');
    tpcall = tpfname(odir,1);
    
    ntp = length(tp);
    C = {};
    for i = 1:ntp
        
        prev = length(C);
        
        % Concatenate transform chain:
        nextname = tp(i).fname;
        C{end+1} = nextname;
        go = true;
        while go
            % Read txt file:
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
            % Determine end of chain:
            if strcmp(nextname,'NoInitialTransform')
                go = false;
            else
                C{end+1} = nextname;
            end
        end
        
        % Write to log file:
        jac = tp(i).jac;
        ncalls = max(jac,length(tp(i).im));
        if ncalls
            if i==1
                permission = 'wt';
            else
                permission = 'at';
            end
            fid = fopen(logfname,permission);
            if fid>2
                fprintf(fid,'\nImages:\n');
                if ~isempty(tp(i).im)
                    fprintf(fid,'   %s\n',tp(i).im{:});
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
        end
        
        % Re-direct previous transform to new level:
        if prev>0
            copytp(tpfname(odir,prev),odir,prev,false);
        end
        
        % Generate copies of all transforms in chain:
        ntot = length(C);
        for j = prev+1:ntot
            copytp(C{j},odir,mod(j-1,ntot)+1,j==ntot);
        end
        pause(0.1);
        
        % Call transformix for selected images/jacobian:
        jac = tp(i).jac;
        for j = 1:max(jac,length(tp(i).im))
            chain = tp(i).chain;
            chlink = sum([tp(1:i).chain]==tp(i).chain);
            if ~isempty(tp(i).im)
                % Fix MHD in case copied:
                fixMHDnames(tp(i).im{j})
                fmt = 'TransformedImage.%02u.%02u.%02u.mhd';
                outfn = fullfile(odir,sprintf(fmt,chain,chlink,j));
                inputs = {'in',tp(i).im{j},'outfn',outfn};
            end
            str = elxobj.sysCmd(odir,inputs{:},'tp',tpcall,'jac',jac,...
                'wait',true,'threads',ncores);
            system(str);
            if jac
                jac = false;
                movefile(fullfile(odir,'spatialJacobian.mhd'),...
                    fullfile(odir,sprintf('spatialJacobian.%02u.%02u.mhd',...
                    chain,chlink)));
                movefile(fullfile(odir,'spatialJacobian.raw'),...
                    fullfile(odir,sprintf('spatialJacobian.%02u.%02u.raw',...
                    chain,chlink)));
            end
            if ~isempty(inputs)
                % Fix copied MHD
                fixMHDnames(outfn);
            end
        end
        
        % Prep for new chain:
        if (i==ntp) || (tp(i+1).chain~=tp(i).chain)
            C = {};
            delete(fullfile(odir,'MultiTransPar.*.txt'));
        end
        
    end
    delete(fullfile(odir,'result.mhd'));
    delete(fullfile(odir,'result.raw'));
    h = toc(t);
end

function str = tpfname(odir,i)
str = fullfile(odir,sprintf('MultiTransPar.%02u.txt',i));

function copytp(fname,odir,i,fin)
% Read file:
    fid = fopen(fname,'r'); str = fread(fid,'*char')'; fclose(fid);
% Change line:
    if fin % End of chain
        itname = 'NoInitialTransform';
    else % Direct to next transform
        itname = tpfname(odir,i+1);
    end
    str = regexprep(str,'\(InitialTransformParametersFileName "([^"]*)',...
        ['\(InitialTransformParametersFileName "',itname]);
% Re-save file:
    fid = fopen(tpfname(odir,i),'w');
    fwrite(fid,str,'char');
    fclose(fid);
    
    