function results = recon4fw(fname1,fname2,id)
% Reconstruct MGE data (MR Solutions or Varian) and save for Fat-Water analysis
% Inputs:
%   flag = T/F for MR Solutions / Varian, respectively (default = Varian)

fprintf('% 2u Inputs\n',nargin);
if nargin==0
    fname1 = pwd;
end
if nargin~=3
    % Select raw MRI data: (Varian)
    fnames = dir(fullfile(fname1,'mge3d*.fid'));
    fnames = {fnames([fnames(:).isdir]).name};
    % Assume mge3d was acquired in pairs: e.g. mge3d_Fat1_01.fid & mge3d_Fat2_01.fid
    nf = floor(length(fnames)/2);
    
    % Ask user for IDs:
    nstr = cellfun(@(x)sprintf('%02u',x),num2cell(1:nf),'UniformOutput',false);
    prompt = strcat(nstr,':');
    defAns = strcat('m',nstr);
    id = inputdlg(prompt,'Set subject ID:',1,defAns);
    if isempty(id)
        return;
    end
    
    fprintf('Starting batch Fat Fraction analysis ...\n');
    for i = 1:nf
        fnames = dir(fullfile(fname1,sprintf('mge3d-Fat*_%02u.fid',i)));
        fnames = fullfile(fname1,{fnames(:).name}');
        fprintf(' ... % 2u : %s\n',i,id{i});
        batch(@recon4fw,1,{fnames{1},fnames{2},id{i}});
        
    end
    fprintf('Done\n');
    results = [];
    return;
elseif ~(ischar(fname1) && ischar(fname2) && ischar(id) ...
        && exist(fname1,'dir') && exist(fname2,'dir'))
    error('Invalid input');
end

fprintf('Starting Fat/Water analysis:\n    %s\n    %s ...\n',fname1,fname2);

imdata = struct('images',[],'TE',[],'FieldStrength',[],'PrecessionIsClockwise',1);
params = struct('species',struct('name',{'water','fat'},...
                                 'frequency',{0,[-3.80, -3.40, -2.60, -1.94, -0.39, 0.60]},...
                                 'relAmps',{1,[0.087 0.693 0.128 0.004 0.039 0.048]}),...
    'size_clique',1,...
    'range_r2star',[0 500],...
    'NUM_R2STARS',51,...
    'range_fm',[-800 800],...
    'NUM_FMS',301,...
    'NUM_ITERS',40,...
    'SUBSAMPLE',2,...
    'DO_OT',1,...
    'LMAP_POWER',2,...
    'lambda',0.05,...
    'LMAP_EXTRA',0.05,...
    'TRY_PERIODIC_RESIDUAL',0);
% params2 = params;
% params2.NUM_MAGN = 1;
% params2.THRESHOLD = 0.04;

% Load kdata:
info = readPROCPAR(fname1);
if ~strcmp(info.seqfil,'mge3d')
    warning('Invalid imaging protocol.')
    return;
end
k = readFID(fname1,info);
info2 = readPROCPAR(fname2);
if ~strcmp(info.seqfil,'mge3d')
    warning('Invalid imaging protocol.')
    return;
elseif (info.tr~=info2.tr) || (info.sw~=info2.sw) || (info.alfa~=info2.alfa) ...
        || (info.B0~=info2.B0)
    warning('Imaging parameters do not match.')
    return;
end
k = cat(4,k,readFID(fname2,info2));
[imdata.TE,ind] = sort([info.TE,info2.TE]);
imdata.FieldStrength = info.B0/10000;
imdata.TE = imdata.TE*1e-3; % needs to be in seconds
k = permute(k(:,:,:,ind,:),[1,2,5,3,4]);
clear info2;

% Reconstruct 3D images:
[d(1),d(2),d(3),d(4),d(5)] = size(k);
imdata.images = zeros(d);
for i = 1:d(5)
    for j = 1:d(4)
        imdata.images(:,:,:,j,i) = fftshift(fftn(k(:,:,:,j,i)));
    end
end

% Calculate Fat/Water Maps using GraphCut:
tdata = imdata;
f = zeros(d(1:3)); w = f; r2s = f; fmap = f;
for i = 1:d(3) % loop over slices since it's a 2D analysis
    % Grab slice:
    tdata.images = imdata.images(:,:,i,:,:);
    
    % Graphcut:
    results = fw_i2cm1i_3pluspoint_hernando_graphcut(tdata,params);
    % Collect results:
    f(:,:,i) = results.species(2).amps;
    w(:,:,i) = results.species(1).amps;
    r2s(:,:,i) = results.r2starmap;
    fmap(:,:,i) = results.fieldmap;
    
    fprintf('Calculating fat slice %u / %u\n',i,d(3));
end
% delete(hw);
results.species(1).amps = w;
results.species(2).amps = f;
results.r2starmap = r2s;
results.fieldmap = fmap;
ff = computeFF(results);
results.fatfraction = ff;

fov = [info.lro,info.lpe,info.lpe2];

outdir = fileparts(fname1);
[outdir,dstr] = fileparts(outdir);
dstr = strsplit(dstr,'_');
dstr = dstr{2};
outdir = fullfile(fileparts(outdir),'Processed',id);
if ~exist(outdir,'dir')
    mkdir(outdir);
end
outname = sprintf('%s_%s',id,dstr);
saveMHD(fullfile(outdir,[outname,'.mhd']),...
    cat(4,abs(w),abs(f),r2s,fmap,ff),...
    {'Water','Fat','R2star','FieldMap','FatFraction'},...
    fov);
save(fullfile(outdir,[outname,'_mge3d.mat']),'imdata','params','results');
