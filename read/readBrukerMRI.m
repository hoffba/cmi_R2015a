function [img,label,fov,p] = readBrukerMRI(varargin)

imname = varargin{1};
[fdir,tname] = fileparts(imname);

% Determine matrix and FOV from parameters:
switch tname
    
    case '2dseq' % Recontructed image
        sdir = fileparts(fileparts(fdir));
        p = readBrukerMRIpar([fullfile(sdir,{'method','acqp'}),{fullfile(fdir,'reco')}]);
        imageObj = ImageDataObject(fdir);
        imageObj.readVisu;
        img = imageObj.data;
        d = size(img);
        fov = p.PVM_Fov; %mm
        
        if ndims(img)>3
            labels4d = cellfun(@(x)strsplit(x,'_'),imageObj.dataDimensions,'UniformOutput',false);
            labels4d = cellfun(@(x)x{end},labels4d,'UniformOutput',false);
        else
            labels4d = {};
        end
        
        clear imageObj;
        
        if strcmp(p.PVM_SpatDimEnum,'<2D>')
            ind = find(strcmp(labels4d,'SLICE'),1);
            ord = 1:length(d);
            ord([3,ind+4]) = [ind+4,3];
            img = permute(img,ord);
            d = d(ord);
            fov(3) = d(3)*p.PVM_SPackArrSliceDistance;
            labels4d(ind) = [];
        end
        img = flip(permute(reshape(img,[d(1:3),prod(d(4:end))]),[2,1,3:length(d)]),1);
        
        if isempty(labels4d)
            label = {p.PULPROG(2:end-5)};
        else
            label = {''};
            for i = 1:length(labels4d)
                switch labels4d{i}
                    case 'DIFFUSION'
                        lstr = strcat('b',cellfun(@num2str,num2cell(p.PVM_DwEffBval),'UniformOutput',false)');
                    case 'ECHO'
                        lstr = strcat('TE',cellfun(@num2str,num2cell(p.EffectiveTE)','UniformOutput',false));
                    otherwise
                        lstr = cellfun(@num2str,num2cell(1:d(i+4)),'UniformOutput',false)';
                end
                nstr = repmat(lstr',length(label),1);
                label = strcat(repmat(label,length(lstr),1),'_',nstr(:));
            end
        end
        
        fov = fov([2,1,3]);
        
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
        fov = fov([2,1,3]);
        
    otherwise
        error('Invalid input file.')
end

% Fix to 4D and adjust labels:
% if d(4)>1
%     ladd = cellfun(@num2str,num2cell(1:d(4)),'UniformOutput',false)';
%     label = strcat(repmat(label,d(4),1),'_rep',ladd);
% end
% if d(5)>1
%     ladd = repmat(cellfun(@num2str,num2cell(1:d(5)),'UniformOutput',false),length(label),1);
%     label = strcat(repmat(label,d(5),1),'_echo',ladd(:));
% end
% if d(6)>1
%     ladd = repmat(cellfun(@num2str,num2cell(1:d(6)),'UniformOutput',false),length(label),1);
%     label = strcat(repmat(label,d(6),1),'_chan',ladd(:));
% end
            
% Image scaling
% for i = 1:nvec
%     if m(i) || b(i)
%         img(:,:,:,i) = img(:,:,:,i) / m(i) - b(i);
%     end
% end
