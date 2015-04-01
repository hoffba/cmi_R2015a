function [kdata,pp] = readFID(fdir,pp)
% read in raw FID data from Varian MRI as k-space
% Input:    fdir = directory containing procpar file
%           pp (optional) = structure containing acquisition parameters
% Output:   kdata = k-space complex matrix
%                   dimensions: [np,nv,ns,ne,narr]
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

        % Load image parameters from "procpar" file
        if strcmp(pp.seqcon{1}(3),'s')
            narr = pp.arraydim/pp.nv;
        else
            narr = pp.arraydim;
        end

        % Reshape kdata into 5D package
        if all(isfield(pp,{'pelist','etl','kzero'}))
            % Determine PE order from table:
            ntrains = pp.nv/pp.etl;
            t1 = reshape(pp.pelist,[pp.etl,ntrains]);
            ind = zeros(pp.etl,pp.ns,ntrains);
            for i = 1:ntrains
                ind(:,:,i) = t1(:,i)*ones(1,pp.ns) + ones(pp.etl,1)*(0:pp.ns-1)*pp.nv;
            end
            [~,ind] = sort(ind(:));
            kdata = reshape(kdata,[ pp.np/2 , pp.nv*pp.ns*pp.ne , pp.arraydim ]);
            kdata = kdata(:,ind,:);

            dims = [pp.np/2 , pp.nv , pp.ns , pp.ne , narr];
            torder = [1,2,3,4,5];
        else
            dims = [pp.np/2 , pp.ne , pp.ns , narr , pp.nv];
            torder = [1,5,3,2,4];
            if strcmp(pp.seqcon{1}(3),'c')
                dims = dims([1,2,3,5,4]);
                torder = [1,4,3,2,5];
            end
        end
        kdata = permute(reshape(kdata,dims),torder);
        
        % Un-interleave if necessary:
        if ~issorted(pp.pss)
            [pp.pss,xi] = sort(pp.pss);
            kdata = kdata(:,:,xi,:,:);
        end
    end
end