function [kdata,fov,label] = cmi_readFID(varargin)
% function to read in k-space complex data from .fid file
% returns kdata in same dimensions as is determined by the acquisition
%     dimensions

kdata = [];
fdir = [];
if (nargin == 0)
    fdir = pwd;
elseif ischar(varargin{1})
    fdir = varargin{1};
end
fname = fullfile(fdir,'fid');

if exist(fname,'file')
    % Read in relevant imaging parameters from "procpar" file
    pars = {'seqfil','lro','lpe','pss','thk','np','nv','ns','ne','arraydim'};
    [seqfil,lro,lpe,pss,thk,np,nv,ns,ne,arraydim]  = cmi_readPROCPAR(fdir,pars);
    if strcmp(seqfil,{'sems'})
        arraydim = arraydim/nv;
    end
    dims = [np/2 ne ns arraydim nv];
    
    % Determine interleaving
    [pss,sind] = sort(pss);
    fov = [lro lpe (pss(end) - pss(1) + thk)];
    
    % Read in complex image data from "fid" file
    fid = fopen(fname,'r','ieee-be');
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
    
    % Reshape kdata into 3D package
    kdata = shiftdim(reshape(kdata,dims),4);
    kdata = permute(kdata(:,:,:,sind,:),[2 1 4 5 3]); % to be [np/2 nv ns arraydim ne]
    
    % Determine labels
    if (arraydim > 1)
        tpar = cmi_readPROCPAR(fdir,{'array'}); % find parameter that was arrayed
        tparvals = cmi_readPROCPAR(fdir,tpar); % find parameter values
        if isempty(tparvals)
            label = regexp(num2str(1:arraydim),'\s*','split');
        else
            label = cell(1,arraydim);
            for i = 1:arraydim
                label{i} = [tpar '=' num2str(tparvals(i))];
            end
        end
    elseif (ne > 1)
        label = regexp(num2str(1:ne),'\s*','split');
    else
        label = {seqfil};
    end
end
