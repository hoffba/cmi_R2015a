% Run Elastix with given image files and list of optimization parameter files
function fnames = elxRegLungsHuman(fnames)

if nargin==0
    fnames.par = strcat('/mnt/cmi/projects/CMI/CMI_Matlab_Programs/BenHoff/cmi/elastix/par0003-lungs/',...
                        {'Par0003.affine.txt','Par0003.bs-R8-ug.txt'});
    [fname,path] = uigetfile('*.mhd','Select Fixed Image:','MultiSelect','off');
    if fname~=0
        fnames.fixed = fullfile(path,fname);
    end
    [fname,path] = uigetfile('*.mhd','Select Moving Image:',path,'MultiSelect','off');
    if fname~=0
        fnames.moving = fullfile(path,fname);
    end
    [fname,path] = uigetfile('*.mhd','Select Fixed VOI:',path,'MultiSelect','off');
    if fname~=0
        fnames.fmask = fullfile(path,fname);
    end
    cd(path);
end

if exist('fnames','var') && all(isfield(fnames,{'fixed','moving','fmask','par'}))
    
    % Check output directory
    if ~isfield(fnames,'out')
        [~,oname,~] = fileparts(fnames.moving);
        fnames.out = fullfile(path,['elxreg_',oname]);
    end
    if ~exist(fnames.out,'dir')
        mkdir(fnames.out);
    end
    
    % Validate file types:
    [~,~,fext] = fileparts(fnames.fixed);
    [~,~,mext] = fileparts(fnames.moving);
    nPar = length(fnames.par);
    pext = cell(1,nPar);
    for i = 1:nPar
        [~,~,pext{i}] = fileparts(fnames.par{i});
    end
    if all(strcmp('.mhd',{fext,mext})) && all(strcmp('.txt',pext))
        % Generate call to Elastix:
        sysstr = ['/opt/elastix/bin/elastix',...
                  ' -f ',fnames.fixed,...
                  ' -m ',fnames.moving,...
                  ' -fMask ',fnames.fmask,...
                  ' -out ',fnames.out,...
                  sprintf(' -p %s',fnames.par{:})];
        system(sysstr);
        % Generate call to Transformix for determinant calculation:
        system(['/opt/elastix/bin/transformix -jac all -out ',fnames.out,...
                ' -tp ',fullfile(fnames.out,'TransformParameters.1.txt')])
    end
end

