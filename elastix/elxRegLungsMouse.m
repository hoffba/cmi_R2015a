% Run Elastix with given image files and list of optimization parameter files
function fnames = elxRegLungsMouse(fnames)

if nargin==0
    [fname,path] = uigetfile('*.mhd','Select Fixed Image:','MultiSelect','off');
    if fname~=0
        fnames(1,:) = {'f',fullfile(path,fname)};
    end
    [fname,path] = uigetfile('*.mhd','Select Moving Image:',path,'MultiSelect','off');
    if fname~=0
        fnames(2,:) = {'m',fullfile(path,fname)};
    end
    [fname,path] = uigetfile('*.mhd','Select Fixed VOI:',path,'MultiSelect','off');
    if fname~=0
        fnames(3,:) = {'fmask',fullfile(path,fname)};
    end
    fnames(4:5,:) = [{'p';'p'},...
                     strcat('/mnt/cmi/projects/CMI/CMI_Matlab_Programs/BenHoff/cmi/elastix/MouseLungs/',...
                            {'MouseLungs.affine.txt','MouseLungs.bs8.txt'})'];
    cd(path);
end

if exist('fnames','var') && iscellstr(fnames)
    
    % Find required file indexes
    indf = find(strcmp('f',fnames(:,1)),1);
    indm = find(strcmp('m',fnames(:,1)),1);
    indfmask = find(strcmp('fmask',fnames(:,1)),1);
    indp = strcmp('p',fnames(:,1));
    
    % Check output directory
    ind = find(strcmp('out',fnames(:,1)),1);
    if isempty(ind)
        fnout = fullfile(fileparts(fnames{indf,2}),'elxreg');
    else
        fnout = fnames{ind,2};
    end
    if ~exist(fnout,'dir')
        mkdir(fnout);
    end
    
    % Validate file types:
    ext = cell(1,size(fnames,1));
    for i = 1:size(fnames,1)
        [~,~,ext{i}] = fileparts(fnames{i,2});
    end
    if ~(isempty(indf) || isempty(indm) || isempty(indfmask)) && any(indp) ...
            && all(strcmp('.mhd',ext(ismember(fnames(:,1),{'f','m','fmask'})))) ...
            && all(strcmp('.txt',ext(ismember(fnames(:,1),{'p'}))))
        try
            % Generate call to Elastix:
            sysstr = ['/opt/elastix/bin/elastix',...
                      ' -f ',fnames{indf,2},...
                      ' -m ',fnames{indm,2},...
                      ' -fMask ',fnames{indfmask,2},...
                      ' -out ',fnout,...
                      sprintf(' -p %s',fnames{indp,2})];
            system(sysstr);
            % Generate call to Transformix for determinant calculation:
%             system(['/opt/elastix/bin/transformix -jac all -out ',fnames.out,...
%                     ' -tp ',fullfile(fnames.out,'TransformParameters.1.txt')])
        catch err
        end
    end
end

