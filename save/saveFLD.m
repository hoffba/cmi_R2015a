function status = saveFLD(fname,img,label,max_ext,info)
%% Save as AVS pseudo-field file
if nargin<5
    info = [];
end
% First check that the file name is correct
[pathstr, name, ext] = fileparts(fname);
if (isempty(ext) || ~strcmp(ext,'.fld'))
    fname = [pathstr name '.fld'];
end
if iscellstr(label)
    for i = 1:length(label)
        label{i}(isspace(label{i})) = '_';
    end
    label = ['label=' sprintf(' "%s"',label{:})];
elseif ischar(label)
    label(isspace(label)) = '_';
    if ~strncmp(label,'label=',6)
        label = ['label= ' label];
    end
else
    label = ['label=' sprintf(' "%d"',1:size(img,4))];
end
img(img<=0) = realmin; % AVS can't handle negative values
min_val = round(double(squeeze(min(min(min(img))))))';
max_val = round(double(squeeze(max(max(max(img))))))';
[avsdim2,avsdim1,ns,veclen] = size(img);
% Slice must be square dimensions
if avsdim1 ~= avsdim2
    tdim = max(avsdim1,avsdim2);
    timg = ones(tdim,tdim,ns,veclen);
    % set padding to min value for each vec
    timg = bsxfun(@times,timg,reshape(min_val,[1 1 1 veclen])); 
    pady = round((tdim - avsdim2)/2);
    padx = round((tdim - avsdim1)/2);
    for i = 1:veclen
        timg((1:avsdim2)+pady,(1:avsdim1)+padx,:,i) = img(:,:,:,i);
    end
    avsdim1 = tdim;
    avsdim2 = tdim;
    max_ext(1:2) = max(max_ext(1:2));
    img = timg;
end

datatyp = 'xdr_short'; % Of output file(s)
fieldtyp = 'rectilinear'; % A purely AVS concept 
min_ext = [0 0 0];
fido = fopen(fname,'w','b');
if fido > 0
    % Check for header info to include:
    if (nargin<5) || ~isfield(info,'RescaleSlope') || isempty(info.RescaleSlope)
        scaleM = num2str(ones(1,veclen));
    end
    if (nargin<5) || ~isfield(info,'RescaleIntercept') || isempty(info.RescaleIntercept)
        scaleB = num2str(zeros(1,veclen));
    end
    if (nargin<5) || ~isfield(info,'AUair') || isempty(info.AUair)
        strAir = [];
    else
        strAir = ['# AUair=' num2str(info.AUair)];
    end
    if (nargin<5) || ~isfield(info,'AUblood') || isempty(info.AUblood)
        strBld = [];
    else
        strBld = ['# AUblood=' num2str(info.AUblood)];
    end
    if (nargin<5) || ~isfield(info,'Modality') || isempty(info.Modality)
        info.Modality = 'Unknown';
    end
    if (nargin<5) || ~isfield(info,'Description') || isempty(info.Description)
        info.Description = '';
    end
    % Write header to AVS field file:
    fprintf(fido,'%s\n',...
        '# AVS field file (@(#)write_field.c     8.1.1.1 Stellar 93/05/12)',...
        ['# TLChenevert matlab-pseudo-AVS file creation date:  ' date],...
        ['# rescaleSlope=' scaleM],...
        ['# rescaleInt=' scaleB]);
    if ~isempty(strAir) && ~isempty(strBld)
        fprintf(fido,'%s\n',...
        strAir,...
        strBld);
    end
    fprintf(fido,'%s\n',...
        ['# modality=' info.Modality],...
        ['# description=' info.Description],...
        'ndim=3	# number of dimensions in the field',...
        ['dim1=' int2str(avsdim1) '	# dimension of axis 1'],...
        ['dim2=' int2str(avsdim2) '	# dimension of axis 2'],...
        ['dim3=' int2str(ns) '	# dimension of axis 3'],...
        'nspace=3	# number of physical coordinates per point',...
        ['veclen=' int2str(veclen) '	# number of components at each point'],...
        ['data=' datatyp '	# portable data format'],...
        ['field=' fieldtyp '	# field type (uniform, rectilinear, irregular)'],...
        label,...
        ['min_ext=' num2str(min_ext) '	# coordinate space extent'],...
        ['max_ext=' num2str(max_ext) '	# coordinate space extent'],...
        ['min_val=' int2str(min_val) '	# min values'],...
        ['max_val=' int2str(max_val) '	# max values'],...
        [12 12]); % required to identify end of header
    status = fseek(fido,-1,'eof');
    if (status == -1)
        input('Error in moving to start of data');
        ferror(fido);
    else
        img = squeeze(permute(img,[4 2 1 3])); % [veclen avsdim2 avsdim1 ns]
        fwrite(fido,img,'short');
        dims = [avsdim1 avsdim2 ns];
        del_ext = (max_ext - min_ext) ./ dims;
        xxx = min_ext(1) + del_ext(1)*((1:avsdim1) - 1);
        yyy = min_ext(2) + del_ext(2)*((1:avsdim2) - 1);
        zzz = min_ext(3) + del_ext(3)*((1:ns) - 1);
        xyz = [xxx yyy zzz];
        status = fseek(fido,0,'eof'); % Confirm we're at end of file
        if (status == -1)
            input('Error in moving to end of data');
            ferror(fido)
        else
            fwrite(fido,xyz,'float');
        end
    end
    status = ~fclose(fido);
else
    status = false;
end