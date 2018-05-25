function [img,label,fov,p] = readBrukerMRI(varargin)

imname = varargin{1};
[fdir,tname] = fileparts(imname);

% Determine matrix and FOV from parameters:
switch tname
    
    case '2dseq' % Recontructed image

        imageObj = ImageDataObject(fdir);
        imageObj.readVisu;
        p = imageObj.Visu;
        img = imageObj.data;
        d = size(img);
        fov = p.VisuCoreExtent;
        
        if ndims(img)>3
            labels4d = cellfun(@(x)strsplit(x,'_'),imageObj.dataDimensions,'UniformOutput',false);
            labels4d = cellfun(@(x)x{end},labels4d,'UniformOutput',false);
        else
            labels4d = {};
        end
        
        clear imageObj;
        
        if p.VisuCoreDim==2
            ind = find(strcmp(labels4d,'SLICE'),1);
            ord = 1:length(d);
            ord([3,ind+4]) = [ind+4,3];
            img = permute(img,ord);
            d = d(ord);
            fov(3) = d(3)*p.VisuCoreSlicePacksSliceDist;
            labels4d(ind) = [];
        end
        img = permute(reshape(img,[d(1:3),prod(d(4:end))]),[2,1,3,4]);
        
        if isempty(labels4d)
            label = {p.VisuAcquisitionProtocol};
        else
            label = {''};
            for i = 1:length(labels4d)
                switch labels4d{i}
                    case 'DIFFUSION'
                        lval = reshape(p.VisuAcqDiffusionBMatrix',3,3,[]);
                        nb = size(lval,3);
                        for j = 1:nb
                            lval(j) = trace(lval(:,:,j));
                        end
                        lval = lval(1:nb);
                        lstr = strcat('b',cellfun(@num2str,num2cell(lval),'UniformOutput',false)');
                    case 'ECHO'
                        lstr = strcat('TE',cellfun(@num2str,num2cell(p.VisuAcqEchoTime)',...
                            'UniformOutput',false));
                    otherwise
                        lstr = cellfun(@num2str,num2cell(1:d(i+4)),'UniformOutput',false)';
                end
                nstr = repmat(lstr',length(label),1);
                label = strcat(repmat(label,length(lstr),1),'_',nstr(:));
            end
        end
        
    case 'fid'  % Raw k-space data
        % * returns img as 5D matrix: [d1,d2,d3,ncoils,d4]
        
        recopart = {'quadrature','phase_rotate','zero_filling','FT'};
        kObj = CKDataObject(fullfile(fdir,'pdata','1'));
        kObj.readReco;
        imObj = kObj.reco(recopart,'image');
        img = flip(permute(imObj.data,[3,1,2,5,6,4]),1);
        fov = imObj.Reco.RECO_fov*10; % mm
        if length(imObj.Method.EffectiveTE)>1
            label = strcat('TE',cellfun(@num2str,num2cell(round(imObj.Method.EffectiveTE*1000)),'UniformOutput',false));
        else
            label = cellfun(@num2str,num2cell(1:size(img,5)),'UniformOutput',false);
        end
        
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
