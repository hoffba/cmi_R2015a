function [img,label,fov] = readMRS(varargin)

if nargin==0
    % Use GUI to select file
    [fname,fpath] = uigetfile('*.SUR;*.MRD','Select a file:');
    [~,~,ext] = fileparts(fname);
else
    [fpath,fname,ext] = fileparts(varargin{1});
    fname = strcat(fname,ext);
end
if ischar(fname)
    % Find files to load:
    allfnames = dir(fullfile(fpath,['*',ext]));
    ind = find(strcmp(fname,{allfnames(:).name}),1);
    [~,dim,par] = Get_mrd_3D1(fullfile(fpath,fname));
    if isfield(par,'no_expts')
        nv = par.no_expts;
    else
        nv = 1;
    end
    if all(isfield(par,{'IM_SLICE','IM_EXPERIMENT','IM_TOTAL_SLICES'}))
        ns = par.IM_TOTAL_SLICES;
        ii = ind - par.IM_SLICE - (par.IM_EXPERIMENT-1)*ns + (1:(ns*nv));
    else
        error('Missing par fields ... fix the code.');
    end
    allfnames = allfnames(ii);
    nf = length(ii);
    img = zeros(dim(1),dim(2),ns,nv);
    hw = waitbar(0,'Reading image slices ...');
    for j = 1:nv
        for i = 1:ns
            ii = i+(j-1)*ns;
            img(:,:,i,j) = Get_mrd_3D1(fullfile(fpath,allfnames(ii).name));
            waitbar(ii/nf,hw);
        end
    end
    delete(hw);
    label = strcat('MRS-',cellfun(@num2str,num2cell(1:nv),'UniformOutput',false));
    fov = par.IM_FOV; fov(3) = fov(3)*ns;
end
