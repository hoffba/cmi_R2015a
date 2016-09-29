function imdata = recon4fw(fname)
% Reconstruct MGE data (MR Solutions or Varian) and save for Fat-Water analysis
% Inputs:
%   flag = T/F for MR Solutions / Varian, respectively (default = Varian)

imdata = struct('images',{[]},'TE',{[]},'FieldStrength',{[]},...
    'PrecessionIsClockwise',{1});
if nargin==0
    [fname,fdir] = uigetfile('*','Select k-space data:');
elseif ~ischar(fname) || ~exist(fname,'file')
    error('Invalid input file name.')
elseif strcmp(fname(end-3:end),'.fid')
    fdir = fname;
    fname = 'fid';
else
    [fdir,fname,ext] = fileparts(fname);
    fname = [fname,ext];
end
if ~ischar(fdir)
    return;
end

switch fname(end-2:end)
    case 'fid'
        [~,fname] = fileparts(fdir);
        info = readPROCPAR(fdir);
        if strcmp(info.seqfil,'mge3d')
            % Read k-space data:
            k = readFID(fdir,info);
            % Organize complex k-space data into 4D format
            k = permute(k,[1,2,5,3,4]);
%             % Assume equally spaced echoes:
%             imdata.TE = (0:(info.ne-1))*info.te2 + info.te;
            imdata.TE = info.TE;
            imdata.FieldStrength = info.B0/10000; % convert Gauss to Tesla
        else
            error('Invalid acquisition protocol. Must be mge3d.');
        end
    case 'MRD'
        [k,info] = Get_mrd_3D1(fullfile(fdir,fname));
        
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % NEED TO FIND THE MRD DATA TO FURTHER IMPLEMENT THIS.
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    otherwise
        error('Invalid file type.');
end

% Reconstruct 3D images:
d = size(k);
imdata.images = zeros(d);
for i = 1:d(5)
    for j = 1:d(4)
        imdata.images(:,:,:,j,i) = fftshift(fftn(k(:,:,:,j,i)));
    end
end

if nargout==0
    % If no output is requested, save to file:
    save(fullfile(fdir,[fname,'_fwData.mat']),'-struct','imdata');
end