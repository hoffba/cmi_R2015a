function varargout = yacta(fnames,varargin)
% Performs YACTA lung processing and analysis
% Syntax: yacta(fname,option1,option2,...)
% Inputs:
%   fname = location of image file for processing
%   options = strings for optional input to yacta64
%       (i.e. 'airways', 'parenchyma', 'renderer', 'hide', 'exportlabels', 'yactascp')
%       ** no optional inputs uses all

if ischar(fnames)
    fnames =  {fnames};
end

nf = numel(fnames);
cmd = cell(nf,1);
for i = 1:nf
    
    fname = fnames{i};
    if isfolder(fname)

        % Compile .dcv file with DICOM names in folder
        [fpath,dcmname] = fileparts(fname);
        tnames = dir(fname);
        tnames([tnames.isdir]) = [];
        tnames = {tnames.name};
        tnames(~cellfun(@isdicom,tnames)) = [];
        fname = fullfile(fname,sprintf('%s.dcv',dcmname));

        % write DCV file
        fid = fopen(fname,'wt');
        fprintf(fid,'%s\n',tnames{:});
        fclose(fid);

    elseif ~ismember(fname(end-3:end),{'.dcv','.mhd'})

        % Must convert to MHD for YACTA to read:
        fpath = fileparts(fname);
        [img,label,fov,orient,info,fnameOut] = cmi_load(1,[],fname);
        fname = fullfile(fpath,[fnameOut,'.mhd']);
        stat = saveMHD(fname,img,label,fov,orient);

    end

    % Validate YACTA options:
    opts = {'airways','parenchyma','renderer','hide','exportlabels','yactascp'};
    if nargin > 1
        opts = varargin(ismember(varargin,opts));
    end

    % Generate system call
    cmd{i} = sprintf('yacta64 "%s"%s',fname,sprintf(' --%s',opts{:}));
    
end

% Send system call:
fprintf('System Call: %s\n',cmd);
t = tic;
stat = system(cmd);
fprintf('Complete (%s)\n',duration(toc(t)));

if nargout
    % Load results to pass back from function:
    varargout = {};
end