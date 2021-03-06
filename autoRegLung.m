function autoRegLung(odir,workingdir,elxfname,f_name,fv_name,m_name,mv_name,varargin)
% Inputs:
%   odir = elastix output directory
%   elxfname = cell or char, Elastix parameter files to use:
%   f_name = fixed image file name (MHD)
%   fv_name = fixed image segmentation file name (MHD)
%   m_name = moving image file name (MHD)
%   mv_name = moving image segmentation file name (MHD)

%% Validate options
    opt = validateInputs(odir,workingdir,elxfname,f_name,fv_name,m_name,mv_name,varargin);

%% Set temp file names:
    [~,tname,~] = fileparts(m_name);
    outfn = [tname,'_R'];

    mo_name = 'elxtemp-origm.mhd';
    ft_name = 'elxtemp-f.mhd';
    mt_name = 'elxtemp-m.mhd';
    fmask_name = 'elxtemp-fMask.mhd';
    t0_name = 'InitialTransform.txt';
    
%% Load images:

    fprintf('Loading fixed image: %s\n',f_name);
    [f_img,~,f_fov,f_info] = readMHD(f_name);
    f_info = struct('dircos',[f_info.SliceOrient,cross(f_info.SliceOrient(1:3),f_info.SliceOrient(4:6))],...
                    'slcpos',f_info.SlicePos);
    f_d = size(f_img);
    f_voxsz = f_fov./f_d;
    f_mask = readMHD(fv_name);
    fv_d = size(f_mask);
    if any(f_d~=fv_d)
        error('Fixed segmentation dimensions do not match image.');
    else
        f_mask = f_mask>0;
    end
    
    fprintf('Loading moving image: %s\n',m_name);
    [m_img,~,m_fov,m_info] = readMHD(m_name);
    m_info = struct('dircos',[m_info.SliceOrient,cross(m_info.SliceOrient(1:3),m_info.SliceOrient(4:6))],...
                    'slcpos',m_info.SlicePos);
    m_d = size(m_img);
    m_voxsz = m_fov./m_d;
    m_mask = readMHD(mv_name);
    mv_d = size(m_mask);
    if any(m_d~=mv_d)
        error('Fixed segmentation dimensions do not match image.');
    else
        m_mask = m_mask>0;
    end
    % Save original moving image for final transform:
    fprintf('Saving original moving image to directory ...\n');
    saveMHD(fullfile(workingdir,mo_name),m_img,{'elxtemp-origm'},m_fov,m_info);

%% Preprocess and save to elastix directory:

    fprintf('Determining initial transform ...\n');
    
    % Automatically find initial transform:
    x = [1 0 0 0 1 0 0 0 1 0 0 0];
    % Scale / Translate based on VOI limits
    [v0y,v0x,v0z] = ind2sub(f_d,find(f_mask));
    [v1y,v1x,v1z] = ind2sub(m_d,find(m_mask));
    v0x = ([min(v0x),max(v0x)] - f_d(1)/2) * f_voxsz(1);
    v0y = ([min(v0y),max(v0y)] - f_d(2)/2) * f_voxsz(2);
    v0z = ([min(v0z),max(v0z)] - f_d(3)/2) * f_voxsz(3);
    v1x = ([min(v1x),max(v1x)] - m_d(1)/2) * m_voxsz(1);
    v1y = ([min(v1y),max(v1y)] - m_d(2)/2) * m_voxsz(2);
    v1z = ([min(v1z),max(v1z)] - m_d(3)/2) * m_voxsz(3);
    x([1,5,9]) = [ diff(v1x)/diff(v0x) ,...
                   diff(v1y)/diff(v0y) ,...
                   diff(v1z)/diff(v0z) ];
    x(isnan(x)) = 1; % correct for divide by zero
    x(10:12) =   [ (sum(v1x)-sum(v0x)*x(1))/2 ,...
                   (sum(v1y)-sum(v0y)*x(5))/2 ,...
                   (sum(v1z)-sum(v0z)*x(9))/2 ] ...
                   + (m_info.slcpos - f_info.slcpos);
    % Save to file:
    defVal = min(0,min(m_img(:)));
    saveTform(fullfile(workingdir,t0_name),x,f_voxsz,f_d,f_info.dircos,f_info.slcpos,defVal);

    % Fixed image:
    fprintf('Preprocessing fixed image ...\n');
    if any(opt.f_filtN)
        if strcmp(opt.f_filtType,'gauss')
            f_img = feval(opt.f_filtHandle,f_img,opt.f_filtN);
        else
            for islc = 1:f_d(3)
                f_img(:,:,islc) = feval(opt.f_filtHandle,f_img(:,:,islc),opt.f_filtN);
            end
        end
    end
    if any(~isinf(opt.f_clim))
        f_img(f_img<opt.f_clim(1)) = opt.f_clim(1);
        f_img(f_img>opt.f_clim(2)) = opt.f_clim(2);
    end
    saveMHD(fullfile(workingdir,ft_name),f_img,{'elxtemp-f'},f_fov,f_info);
    clear f_img
    
    % Moving image:
    fprintf('Preprocessing moving image ...\n');
    if any(opt.m_filtN)
        if strcmp(opt.m_filtType,'gauss')
            m_img = feval(opt.m_filtHandle,m_img,opt.m_filtN);
        else
            for islc = 1:f_d(3)
                m_img(:,:,islc) = feval(opt.m_filtHandle,m_img(:,:,islc),opt.m_filtN);
            end
        end
    end
    if any(~isinf(opt.f_clim))
        m_img(m_img<opt.m_clim(1)) = opt.m_clim(1);
        m_img(m_img>opt.m_clim(2)) = opt.m_clim(2);
    end
    saveMHD(fullfile(workingdir,mt_name),m_img,{'elxtemp-m'},m_fov,m_info);
    clear m_img

    % Fixed mask:
    fprintf('Preprocessing fixed mask ...\n');
    if any(opt.f_dilateN)
        r = opt.f_dilateN;
        ni = max(ceil(r/5)); % 5 voxel radius per step max, otherwise takes way too long
        for j = 1:ni
            rt = min(r,5);
            r = r - rt;
            se = bwellipsoid(rt);
            f_mask = imdilate(f_mask,se);
        end
    end
    saveMHD(fullfile(workingdir,fmask_name),f_mask,{'elxtemp-fmask'},f_fov,f_info);
    clear f_mask


%% Call Elastix/Transformix:
% Copy Elastix parameter files to output directory:
for i = 1:length(elxfname)
    copyfile(elxfname{i},workingdir);
    [~,elxfname{i}] = fileparts(elxfname{i});
end
% Generate system call:
optstr = '';
if opt.jac
    optstr = [optstr,' -jac all'];
end
if opt.jacmat
    optstr = [optstr,' -jacmat all'];
end
if opt.def
    optstr = [optstr,' -def all'];
end
finalT_name = sprintf('TransformParameters.%u.txt',length(elxfname)-1);
txferstr = '';
if ~strcmp(odir,workingdir)
    txferstr = sprintf(' && cp -a %s/. %s',workingdir,odir);
end

% Set environment variables:
str = sprintf('export odir="%s"',workingdir);
% Call to Elastix:
str = [str,sprintf(' && elastix -f $odir"/%s" -m $odir"/%s" -out $odir -fMask $odir"/%s" -t0 $odir"/%s"',...
                   ft_name,mt_name,fmask_name,t0_name)];
str = [str,sprintf(' -p $odir"/%s.txt"',elxfname{:})];
% Call to Transformix for final transform:
str = [str,sprintf(' && transformix -out $odir -tp $odir"/%s" -in $odir"/%s" ',finalT_name,mo_name),optstr];
% Rename transformed result:
str = [str,sprintf(' cp $odir"/result.mhd" $odir"/%1$s.mhd" && cp $odir"/result.raw" $odir"/%1$s.raw"',outfn)];
% Delete temporary files:
str = [str,' && find $odir -name "elxtemp-*" -exec rm -f {} \;'];
str = [str,' && find $odir -name "result.*" -exec rm -f {} \;'];
% Transfer back to Turbo:
str = [str,txferstr];
% Delete temporary scratch folder: 
str = [str,' && rm -r $odir;'];

% Save full command string for debug:
fid = fopen(fullfile(workingdir,'elastixCMD.txt'),'w');
if fid
    fprintf(fid,'%s',str);
    fclose(fid);
end

fprintf('Sending Elastix call ...\n');
system(str);

fprintf('DONE!\n');


function opt = validateInputs(odir,workingdir,elxfname,f_name,fv_name,m_name,mv_name,opts)
% Validate input parameters
p = inputParser;

addRequired(p,'odir',@ischar);
addRequired(p,'workingdir',@ischar);
addRequired(p,'elxfname',@(x)iscellstr(x));
addRequired(p,'f_name',@(x)exist(x,'file'));
addRequired(p,'fv_name',@(x)exist(x,'file'));
addRequired(p,'m_name',@(x)exist(x,'file'));
addRequired(p,'mv_name',@(x)exist(x,'file'));

addParameter(p,'f_clim',[-inf,inf],@(x)isnumeric(x)&&numel(x)==2&&x(1)<x(2));
addParameter(p,'f_filtType','wiener',@(x)ischar(x)&&ismember(x,{'gauss','wiener','median'}));
addParameter(p,'f_filtN',zeros(1,3),@(x)isnumeric(x)&&numel(x)==3);

addParameter(p,'f_dilateN',zeros(1,3),@(x)isnumeric(x)&&numel(x)==3);

addParameter(p,'m_clim',[-inf,inf],@(x)isnumeric(x)&&numel(x)==2&&x(1)<x(2));
addParameter(p,'m_filtType','wiener',@(x)ischar(x)&&ismember(x,{'gauss','wiener','median'}));
addParameter(p,'m_filtN',zeros(1,3),@(x)isnumeric(x)&&numel(x)==3);

addParameter(p,'jac',false,@(x)islogical(x)&&numel(x)==1);
addParameter(p,'jacmat',false,@(x)islogical(x)&&numel(x)==1);
addParameter(p,'def',false,@(x)islogical(x)&&numel(x)==1);
addParameter(p,'Nthreads',4,@(x)isnumeric(x)&&numel(x)==1);

parse(p,odir,workingdir,elxfname,f_name,fv_name,m_name,mv_name,opts{:});
opt = p.Results;

opt.f_filtN = round(opt.f_filtN);
opt.f_dilateN = round(opt.f_dilateN);
opt.m_filtN = round(opt.m_filtN);
opt.Nthreads = round(opt.Nthreads);

opt.f_filtN = 2*opt.f_filtN + 1;
switch opt.f_filtType
    case 'gauss'
        opt.f_filtN = opt.f_filtN/3;
        opt.f_filtHandle = @(x,f)imfilter(x,gaussND(f));
    case 'wiener'
        opt.f_filtHandle = @(x,f)wiener2(x,f(1:2));
    case 'median'
        opt.f_filtHandle = @(x,f)medfilt2(x,f(1:2));
end
opt.m_filtN = 2*opt.m_filtN + 1;
switch opt.m_filtType
    case 'gauss'
        opt.m_filtN = opt.m_filtN/3;
        opt.m_filtHandle = @(x,f)imfilter(x,gaussND(f));
    case 'wiener'
        opt.m_filtHandle = @(x,f)wiener2(x,f(1:2));
    case 'median'
        opt.m_filtHandle = @(x,f)medfilt2(x,f(1:2));
end

if ~isdir(opt.odir)
    mkdir(opt.odir);
end
if ~isdir(opt.workingdir)
    mkdir(opt.workingdir);
end

function saveTform(fname,x,voxsz,d,dircos,pos,defVal)
% Saves initial transform parameter file:
t = struct(...
           'NumberOfParameters',length(x),...
           'TransformParameters',x,...
           'Transform','AffineTransform',...
           'InitialTransformParametersFileName','NoInitialTransform',...
           'HowToCombineTransforms','Compose',...
           'FixedImageDimension',3,...
           'MovingImageDimension',3,...
           'Size',d,...
           'Index',zeros(1,3),...
           'Spacing',voxsz,...
           'Origin',pos,...
           'CenterOfRotationPoint',zeros(1,3),...
           'Direction',dircos,...
           'UseDirectionCosines','true',...
           'ResampleInterpolator','FinalBSplineInterpolator',...
           'FinalBSplineInterpolationOrder',3,...
           'Resampler','DefaultResampler',...
           'DefaultPixelValue',defVal,...
           'ResultImageFormat','mhd',...
           'ResultImagePixelType','float',...
           'FixedInternalImagePixelType','float',...
           'MovingInternalImagePixelType','float',...
           'CompressResultImage','false');

fid = fopen(fname,'w');
if fid>2
    fstr = fieldnames(t);
    for i = 1:length(fstr)
        val = t.(fstr{i});
        vstr = '';
        if ischar(val)
            vstr = [' "',val,'"'];
        elseif iscellstr(val)
            vstr = sprintf(' "%s"',val{:});
        elseif all(round(val)==val)
            vstr = sprintf(' % .0f',val);
        elseif isnumeric(val)
            vstr = sprintf(' % .8f',val);
        end
        fprintf(fid,'(%s%s)\n',fstr{i},vstr);
    end
    fclose(fid);
else
    error(['Could not open file: ',fname])
end

