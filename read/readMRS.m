function [img,label,fov] = readMRS(varargin)
% .MRD = Raw k-space data
% .SUR = Reconstructed images

if nargin==0
    % Use GUI to select file
    [fname,fpath] = uigetfile('*.SUR;*.MRD','Select a file:');
    [~,~,ext] = fileparts(fname);
else
    [fpath,fname,ext] = fileparts(varargin{1});
    fname = strcat(fname,ext);
end
if ischar(fname)
    [~,label,~] = fileparts(fname);
    label = strtok(label,'_');
    if strcmpi(ext,'.SUR') % Magnitude images
        
        % Find all files to load:
        fname = dir(fullfile(fpath,[label,'_*',ext]));
        nf = length(fname);
        
        % Load slices:
        [timg,par] = Get_mrd_3D1(fullfile(fpath,fname(1).name));
        img = nan([size(timg),nf]);
        img(:,:,1) = timg;
        pos = nan(nf,3);
        pos(1,:) = [par.IM_EXPERIMENT,par.IM_ECHO,par.IM_SLICE];
        for i = 2:nf
            [img(:,:,i),par] = Get_mrd_3D1(fullfile(fpath,fname(i).name));
            pos(i,:) = [par.IM_EXPERIMENT,par.IM_ECHO,par.IM_SLICE];
        end
        
        % Rearrange and sort:
        [~,ind] = sortrows(pos);
        img = reshape(img(:,:,ind),par.IM_RESOLUTION(1),...
            par.IM_RESOLUTION(2),par.IM_TOTAL_SLICES,[]);
        
        fov = par.IM_FOV;
        
    elseif strcmpi(ext,'.MRD') % Complex k-space data
        
        % Read all k-space data
        [img,par] = Get_mrd_3D1(fullfile(fpath,fname));
        
        % Rearrange into 4D:
        [d(1),d(2),d(3),d(4)] = size(img);
        img = reshape(img,d);
        
        fov = par.FOV;
    else
        error('Invalid file name.');
    end
    
    % Determine FOV:
    if isfield(par,'no_views_2') && (par.no_views_2>1) % 3D acquisition
        fovz = par.SLICE_THICKNESS;
    else
        fovz = max(par.FOV_OFFSETS(:,3)) - min(par.FOV_OFFSETS(:,3))...
            + diff(par.FOV_OFFSETS(1:2,3));
    end
    fov = [ fov , fov , fovz ];
    
    % Determine labels:
    d4 = size(img,4);
    label = repmat({label},1,size(img,4));
    if isfield(par,'no_echoes') && (par.no_echoes>1)
        val = repmat((1:par.no_echoes)',1,d4/par.no_echoes);
        val = cellfun(@(x)num2str(x,'%02u'),num2cell(val(:)),'UniformOutput',false)';
        label = strcat(label,'_Echo',val);
    end
    if isfield(par,'no_experiments') && (par.no_experiments>1)
        val = repmat((1:par.no_echoes),d4/par.no_experiments,1);
        val = cellfun(@(x)num2str(x,'%02u'),num2cell(val(:)),'UniformOutput',false)';
        label = strcat(label,'_Exp',val);
    end
    
end
