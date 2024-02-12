function varargout = yacta(fnames,varargin)
% Performs YACTA lung processing and analysis
% Syntax: yacta(fname,option1,option2,...)
% Inputs:
%   fname = location of image file for processing
%   options = strings for optional input to yacta64
%       ('airways', 'parenchyma', 'renderer', 'hide', 'exportlabels', 'yactascp')
%           ** no optional inputs uses all
%       ('wait' option to suspend action until completion of system call)
% ** Options for YACTA are set in the yacta.ini file
%       (typically found in C:\Users\<*username*>\AppData\Roaming\yacta64)

if ischar(fnames)
    fnames =  {fnames};
end

waitstr = '&';
ind = strcmp('wait',varargin);
if nnz(ind)
    waitstr = '';
    varargin(ind) = [];
end

% Validate YACTA options:
opts = {'airways','parenchyma','renderer','hide','exportlabels','yactascp'};
if numel(varargin)
    opts = varargin(ismember(varargin,opts));
end
if numel(opts)
    optstr = sprintf(' --%s',opts{:});
else
    optstr = '';
end

nf = numel(fnames);
cmd = cell(nf,1);
for i = 1:nf
    
    fname = fnames{i};
    if isfolder(fname)

        % Compile .dcv file with DICOM names in folder
        [~,dcmname] = fileparts(fname);
        tnames = dir(fname);
        tnames([tnames.isdir]) = [];
        tnames = {tnames.name};
        tnames(~cellfun(@(x)exist(x,'file')&&isdicom(x),tnames)) = [];
        svname = fullfile(fname,sprintf('%s.dcv',dcmname));

        % write DCV file
        fid = fopen(svname,'wt');
        fprintf(fid,'%s\n',tnames{:});
        fclose(fid);

    elseif ~ismember(fname(end-3:end),{'.dcv','.mhd'})

        % Must convert to MHD for YACTA to read:
        [img,label,fov,orient,~,~] = cmi_load(1,[],fname);
        [fpath,fname,ext] = fileparts(fname);
        if strcmp(ext,'.gz')
            [~,fname] = fileparts(fname);
        end
        fname = fullfile(fpath,[fname,'.mhd']);
        stat = saveMHD(fname,img,label,fov,orient);

    end

    % Generate system call
    yactastr = 'C:\Program Files\yacta64\yacta64.exe';
    cmd{i} = sprintf('"%s" "%s"%s',yactastr,fname,optstr);
    
end

if nf > 1
    % Run .bat file:
    batname = fullfile(tempdir,sprintf('yacta_%s.bat',datestr(datetime,'yyyymmddHHMMss')));
    fid = fopen(batname,'wt');
    fprintf(fid,'%s\n',cmd{:});
    fclose(fid);
    cmd = batname;
else
    cmd = cmd{1};
end

% Send system call:
fprintf('System Call: %s\n',cmd);
t = tic;
stat = system(sprintf('%s%s',cmd,waitstr));
if isempty(waitstr)
    fprintf('Complete (%s)\n',duration(0,0,toc(t)));
end

if nargout
    % Load results to pass back from function:
    varargout = stat;
end