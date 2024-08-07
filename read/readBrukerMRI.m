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

        
        
        
        
        
        [~,id] = fileparts(fdir);
        p = readBrukerMRIpar(fullfile(fdir,{'acqp','method'}));
        fov = p.PVM_Fov; % mm
        
        switch p.AQ_mod
            case ('qf')
                isComplex = false;
            case ('qseq')
                isComplex = true;
            case ('qsim')
                isComplex = true;
            case ('qdig')
                isComplex = true;
            otherwise
                error('The value of parameter AQ_mod is not supported');
        end
        
        % [ RO x Nphase x Nslices x Nreps x Nechoes x Ncoils ]
        ncoils = p.PVM_EncNReceivers;
        d = p.PVM_Matrix;
        a = p.ACQ_size;
        nechoes = p.PVM_NEchoImages;
        nrep = p.PVM_NRepetitions;
        d(4) = nechoes * nrep * ncoils;
        if d(3)==0
            chk3d = false;
            d(3) = p.PVM_SPackArrNSlices;
            a(3) = d(3);
            fov(3) = (p.PVM_SliceThick+p.PVM_SPackArrSliceGap)*d(3);
        else
            chk3d = true;
        end

        switch p.GO_raw_data_format
            case 'GO_32BIT_SGN_INT'
                f = 'int32';
                bits = 32;
            case 'GO_16BIT_SGN_INT'
                f = 'int16';
                bits = 16;
            case 'GO_32BIT_FLOAT'
                f = 'float32';
                bits = 32;
            otherwise
                f = 'int32';
                bits = 32;
        end
        
        switch p.BYTORDA
            case 'little'
                bo = 'l';
            case 'big'
                bo = 'b';
            otherwise
                bo = 'n';
        end
        
        switch p.GO_block_size
            case 'continuous'
                sz = bits * prod(d) * (1+isComplex);
            case 'Standard_KBlock_Format'
                fact = bits/8 / 1024; % padded to 1kB blocks
                blksz = ceil(a(1)*ncoils*fact)/fact;
                sz = [ blksz , prod(a(2:end))*nrep*nechoes ];
            otherwise
        end
        
        % Read binary data file:
        fid = fopen(imname,'r');
        img = fread(fid,sz,f,0,bo);
        fclose(fid);
        
        % Remove zero-lines:
        if blksz ~= a(1)*ncoils
            img = reshape(img(1:(a(1)*ncoils),:),[]);
        else
        end
        
        if isComplex
            img = complex(img(1:2:end,:),img(2:2:end,:));
            a(1) = a(1)/2;
        end
        
        % Adjust to image matrix:
        img = reshape(permute(reshape(img,[a(1),ncoils,nechoes,a(2:3),nrep]),[1,4,5,2,3,6]),[a,d(4)]);
        padsz = zeros(1,4);
        % Partial Fourier on PE
        [~,ix] = sort(p.ACQ_spatial_phase_1);
        img = img(:,ix,:,:);
        padsz(2) = d(2) - a(2);
        % Partial Fourier on PE2
        [~,ix] = sort(p.ACQ_spatial_phase_2);
        img = img(:,:,ix,:);
        padsz(3) = d(3) - a(3);
        if any(padsz>0)
            img = padarray(img,padsz,0,'pre');
        end
        
        if chk3d
            for i = 1:d(4)
                img(:,:,:,i) = fftshift(fftn(img(:,:,:,i)));
            end
        else
            for i = 1:d(3)
                for j = 1:d(4)
                    img(:,:,i,j) = fftshift(fft2(img(:,:,i,j)));
                end
            end
            if strcmp(p.PVM_ObjOrderScheme,'Interlaced')
                img = img(:,:,p.PVM_ObjOrderList+1,:);
            end
        end

        label = strcat(id,'_',p.ACQ_protocol_name(2:end-1),'_',...
            cellfun(@num2str,num2cell(1:d(4)),'UniformOutput',false));
        
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
