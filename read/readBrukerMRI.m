function [img,label,fov,p] = readBrukerMRI(varargin)

imname = varargin{1};
[fdir,tname] = fileparts(imname);
img = [];
label = {};
fov = [];

% Determine matrix and FOV from parameters:
switch tname
    
    case '2dseq' % Recontructed image
        sdir = fileparts(fileparts(fdir));
        pfnames = [fullfile(sdir,{'method','acqp'}),{fullfile(fdir,'reco')}];
        pfnames(~cellfun(@exist,pfnames)) = [];
        p = readBrukerMRIpar(pfnames);
        imageObj = ImageDataObject(fdir);
        imageObj.readVisu;
        img = imageObj.data;
        [d(1),d(2),d(3),d(4),d(5),d(6)] = size(img);
        fov = p.PVM_Fov; %mm
        
        if isempty(imageObj.FrameGroupParams)
            labels4d = {imageObj.Visu.VisuAcquisitionProtocol};
        elseif isfield(imageObj.FrameGroupParams,'VisuAcqDiffusionBMatrix')
            ldim = cellfun(@(x)~isempty(strfind(x,'DIFFUSION')),imageObj.dataDimensions);
            ddim = find(ldim,1);
            ind = num2cell(ones(1,length(ldim)));
            ind{ddim} = 1:size(imageObj.FrameGroupParams,ddim);
            labels4d = [imageObj.FrameGroupParams(ind{:}).VisuFGElemComment];
        elseif isfield(imageObj.FrameGroupParams,'VisuFGElemId')
            labels4d = [imageObj.FrameGroupParams(:,1).VisuFGElemId];
        else
            labels4d = cellfun(@num2str,num2cell(1:d(4)),'UniformOutput',false);
        end
        
        if strcmp(p.PVM_SpatDimEnum,'<2D>')
            dimstr = cellfun(@(x)strsplit(x,'_'),imageObj.dataDimensions,'UniformOutput',false);
            dimstr = cellfun(@(x)x{end},dimstr,'UniformOutput',false);
            ind = find(strcmp(dimstr,'SLICE'),1);
            ord = 1:6;
            ord([3,ind+4]) = [ind+4,3];
            img = permute(img,ord);
            d = d(ord);
            fov(3) = d(3)*p.PVM_SPackArrSliceDistance;
        end
        img = flip(permute(reshape(img,[d(1:3),prod(d(4:end))]),[2,1,3:length(d)]),1);
        
        if isempty(labels4d)
            label = {p.PULPROG(2:end-5)};
        else
            label = labels4d;
        end
        
    case 'fid'  % Raw k-space data
        % * returns img as 5D matrix: [d1,d2,d3,ncoils,d4]
        
        p = readBrukerMRIpar(fullfile(fdir,{'method','acqp'}));
        recopart = {'quadrature','phase_rotate','zero_filling','FT'};
        kObj = CKDataObject(fullfile(fdir,'pdata','1'));
        kObj.readReco;
        imageObj = kObj.reco(recopart,'image');
        img = flip(permute(imageObj.data,[1,2,3,5,6,4]),1);
        fov = imageObj.Reco.RECO_fov*10; % mm
        if length(imageObj.Method.EffectiveTE)>1
            label = strcat('TE',cellfun(@num2str,num2cell(round(imageObj.Method.EffectiveTE*1000)),'UniformOutput',false));
        else
            label = cellfun(@num2str,num2cell(1:size(img,5)),'UniformOutput',false);
        end
        d = size(img);
        
    otherwise
        error('Invalid input file.')
end

% Bruker position listed as center of image, need change to corner:
p.SlicePos = [p.PVM_SPackArrPhase1Offset,p.PVM_SPackArrReadOffset,p.PVM_SPackArrSliceOffset] - fov.*(1-1./d(1:3))/2;
% Re-orient to xyz:
sporient = p.PVM_SPackArrGradOrient;
dorder = ones(1,3);
for i = 1:3
    dorder(i) = find(sporient(:,:,i)==max(sporient(:,:,i)));
end
[~,dorder] = sort(dorder);
fov = fov(dorder);
p.SlicePos = p.SlicePos(dorder);
p.SliceOrient = reshape(p.PVM_SPackArrGradOrient(:,:,dorder),1,[]);
img = permute(img,[dorder,4:ndims(img)]);

