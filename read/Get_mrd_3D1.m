% Description: Function to open multidimensional mrd images given a filename with PPR-parsing
% Inputs: string filename, double xsize, double ysize
% Outputs: complex data, raw dimension [no_expts,no_echoes,no_slices,no_views,no_views_2,no_samples], MRD/PPR parameters
% Author: Ruslan Garipov
% Date: 01/03/2014 - swapped views and views2 dimension - now correct

function [img,dim,par] = Get_mrd_3D1(filename,reorder1,reorder2)

if nargin<2
    reorder1 = false;
end
if nargin<3
    reorder2 = false;
end


% Read in wrapped imaged and sensitivity maps
fid = fopen(filename,'r');      % Define the file id
val = fread(fid,4,'int32');
xdim = val(1);
ydim = val(2); 
zdim = val(3);
dim4 = val(4);
fseek(fid,18,'bof');
datatype=fread(fid,1, 'uint16');
datatype = dec2hex(datatype);
fseek(fid,48,'bof');
% scaling=fread(fid,1, 'float32');
% bitsperpixel=fread(fid,1, 'uchar');
fseek(fid,152,'bof');
val = fread(fid,2, 'int32');
dim5 = val(1);
dim6 = val(2);
fseek(fid,256,'bof');
% text=fread(fid,256);
no_samples = xdim;
no_views = ydim;
no_views_2=zdim;
no_slices = dim4;
no_echoes = dim5;
no_expts = dim6;

if size(datatype,2)>1
    onlydatatype = datatype(2);
    iscomplex = true;
else
    onlydatatype = datatype(1);
    iscomplex = false;
end
switch onlydatatype
    case '0' 
        dataformat = 'uchar';   %datasize = 1; % size in bytes
    case '1' 
        dataformat = 'schar';   %datasize = 1; % size in bytes
    case '2' 
        dataformat = 'short';   %datasize = 2; % size in bytes
    case '3' 
        dataformat = 'int16';   %datasize = 2; % size in bytes
    case '4' 
        dataformat = 'int32';   %datasize = 4; % size in bytes
    case '5' 
        dataformat = 'float32'; %datasize = 4; % size in bytes
    case '6' 
        dataformat = 'double';  %datasize = 8; % size in bytes
    otherwise
        dataformat = 'int32';   %datasize = 4; % size in bytes
end

dim = [no_expts , no_echoes , no_slices , no_views , no_views_2 , no_samples];
ntot = prod(dim)*(iscomplex+1);
fseek(fid,512,'bof');
[img,count] = fread(fid,ntot,dataformat); % reading all the data at once

if (count~=ntot)
    error('Incomplete data file ...');
end

if iscomplex
    img = complex(img(1:2:end),img(2:2:end));
end

dorder = [6,4,5,3,2,1];
img = permute(reshape(img,dim),dorder);
dim = dim(dorder);
if reorder1
    img = img(:,reshape([no_views/2+1:no_views ; 1:no_views/2],1,no_views),:,:,:,:);
end
if reorder2
    img = img(:,:,reshape([no_views2/2+1:no_views2 ; 1:no_views2/2],1,no_views2),:,:,:);
end
% img = permute(reshape(img,dim),[1,2,3,5,4,6]);
% dim = dim([1,2,3,5,4,6]);
% if reorder1
%     img = img(:,:,:,:,reshape([no_views/2+1:no_views ; 1:no_views/2],1,no_views),:);
% end
% if reorder2
%     img = img(:,:,:,reshape([no_views2/2+1:no_views2 ; 1:no_views2/2],1,no_views2),:,:);
% end
% img = permute(img,6:-1:1);

fseek(fid,120,0);
ppr_text = char(fread(fid,Inf,'uchar')');
fclose(fid);

% parse fields in ppr section of the MRD file
ppr_text = strsplit(ppr_text(2:end),[char([13,10]),':'])';
par.filename = filename;
if numel(ppr_text)>0
    for i = 1:length(ppr_text)
        [str,rem] = strtok(ppr_text{i},' ');
        if ~isempty(rem)
            rem(1) = [];
            switch str
                case {'DISCARD','EXPERIMENT_ARRAY','FOV_READ_OFF','FOV_PHASE_OFF',...
                      'FOV_SLICE_OFF','NO_AVERAGES','NO_ECHOES','NO_RECEIVERS',...
                      'NO_SAMPLES','NO_SLICES','NO_VIEWS','NO_VIEWS_2','PHASE_CYCLE',...
                      'SLICE_BLOCK','SLICE_INTERLEAVE','VIEWS_PER_SEGMENT',...
                      'PHASE_ORIENTATION','X_ANGLE','Y_ANGLE','Z_ANGLE','VAR',...
                      'SAMPLE_PERIOD','SAMPLE_PERIOD_2'}
                    str = strsplit(rem,', ');
                    par.(str{1}) = str2double(str{2});
                case {'IM_FOV','IM_FFT','IM_ORIENTATION','IM_FREQDIR',...
                      'IM_OFFSETS','IM_SLICE','IM_TOTAL_SLICES','IM_ECHO',...
                      'IM_EXPERIMENT','IM_RESOLUTION'}
                    par.(str) = cellfun(@str2double,strsplit(rem,' '));
                case {'OBSERVE_FREQUENCY','PPL'}
                    par.(str) = rem;
                case {'FOV','SMX','SMY','SWX','SWY','SMZ','SWZ'}
                    par.(str) = str2double(rem);
                case {'SLICE_THICKNESS','SLICE_SEPARATION','VAR_ARRAY'}
                    str = strsplit(rem,', ');
                    par.(str{1}) = cellfun(@str2double,str(2:end));
            end
        end
    end
    if isfield('OBSERVE_FREQUENCY','par')
        par.Nucleus = strtok(par.OBSERVE_FREQUENCY,'"');
    else
        par.Nucleus = 'Unspecified';
    end
    par.datatype = datatype;
    file_pars = dir(filename);
    par.date = file_pars.date;
end
