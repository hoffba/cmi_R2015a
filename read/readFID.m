function [kdata,pp] = readFID(fdir,pp)
% read in raw FID data from Varian MRI as k-space
% Input:    fdir = directory containing procpar file
%           pp (optional) = structure containing acquisition parameters
% Output:   kdata = k-space complex matrix
%                   dimensions: [np,nv,ns,narray]
%           pp = structure containing acquisition parameters

kdata = [];
if ischar(fdir) && exist(fullfile(fdir,'fid'),'file')
    if nargin<2
        pp = readPROCPAR(fdir);
    end
        
    % Read in complex image data from "fid" file
    fid = fopen(fullfile(fdir,'fid'),'r','ieee-be');
    if fid>0
        nblocks   = fread(fid,1,'int32');
        ntraces   = fread(fid,1,'int32'); %
        np        = fread(fid,1,'int32'); % readout*2
        fseek(fid,14,0);
        status    = fread(fid,1,'int16');
        fseek(fid,4,0);
        if bitget(status,3)
            dtype = 'int32';
        elseif bitget(status,4)
            dtype = 'float32';
        else
            dtype = 'int16';
        end
        kdata = fread(fid,[np*ntraces+7,nblocks],dtype);
        fclose(fid);
        kdata = complex(kdata(8:2:end,:),kdata(9:2:end,:));

        % Determine data order:
        % seqcon: (1) = multi-echo
        %         (2) = slices
        %         (3) = PE
        %         (4) = PE2
        %         (5) = PE3
        seqcon = pp.seqcon{1};
        np = pp.np/2; % # of complex RO points
        ne = pp.ne; % # of echoes (always compressed)
        nv = pp.nv; % PE steps
        nv2 = pp.nv2;
        nv3 = pp.nv3;
        ns = max(pp.ns,length(pp.pss)); % # of slices
        narr = pp.arraydim; % Arrayed dimension
        arrstr = strsplit(pp.array{1},',');
        ntrace = [np ne];
        nblock = narr;
        pdimt = [1,6];
        pdimb = 7;
        arri = 1;
        if strcmp(seqcon(2),'s')
            nblock = [ns,nblock];
            pdimb = [5,pdimb];
            arrstr(end) = []; % Remove pss from array list
            arri = 2;
            nblock(arri) = nblock(arri)/ns;
        else
            ntrace = [ntrace,ns];
            pdimt = [pdimt,5];
        end
        if strcmp(seqcon(3),'s')
            nblock = [nblock,nv];
            nblock(arri) = nblock(arri)/nv;
            pdimb = [pdimb,2];
        else
            ntrace = [ntrace,nv];
            pdimt = [pdimt,2];
        end
        if strcmp(seqcon(4),'s')
            nblock = [nblock,nv2];
            pdimb = [pdimb,3];
            nblock(arri) = nblock(arri)/nv2;
        else
            ntrace = [ntrace,nv2];
            pdimt = [pdimt,3];
        end
        if strcmp(seqcon(5),'s')
            nblock = [nblock,nv3];
            pdimb = [pdimb,4];
            nblock(arri) = nblock(arri)/nv3;
        else
            ntrace = [ntrace,nv3];
            pdimt = [pdimt,4];
        end
        readdim = [ntrace,nblock];
        readdim = max(readdim,1);
        [~,pdim] = sort([pdimt,pdimb]);
        kdata = permute(reshape(kdata,readdim),pdim);
        d = readdim(pdim);
        % Order should now be:
        %   [np nv nv2 nv3 ns ne narr]

        % Reshape kdata into 5D package
        if all(isfield(pp,{'pelist','etl','kzero'}))
            % Determine PE order from table:
            ntrains = pp.nv/pp.etl;
            t1 = reshape(pp.pelist,[pp.etl,ntrains]);
            ind = zeros(pp.etl,ntrains);
            for i = 1:ntrains
                ind(:,:,i) = t1(:,i) + ones(pp.etl,1)*pp.nv;
            end
            [~,ind] = sort(ind(:));
            kdata = kdata(:,ind,:,:,:,:,:);
        end
        
        % Un-interleave if necessary:
        if ~issorted(pp.pss)
            [pp.pss,ind] = sort(pp.pss);
            kdata = kdata(:,:,:,:,ind,:,:);
        end
        
        % Determine final dimensions in 4D:
        lbl4D = {};
        if nv2==0
            % put slices in 3rd dimension
            kdata = permute(kdata,[1,2,5,6,7,3,4]);
        else % 3D imaging
            if ns>1
                lbl4D = [lbl4D,{'Slice'}];
            end
            % probably never need nv3, so put at end
            kdata = permute(kdata,[1,2,3,5,6,7,4]);
        end
        if ne>1
            lbl4D = [lbl4D,{'Echo'}];
        end
        if d(7)>1
            lbl4D = [lbl4D,arrstr];
        end
        pp.array = lbl4D;
        [d(1),d(2),d(3),~] = size(kdata);
        kdata = reshape(kdata,d(1),d(2),d(3),[]);
    end
end