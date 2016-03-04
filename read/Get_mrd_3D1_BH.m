% Description: Function to open multidimensional mrd images given a filename with PPR-parsing
% Inputs: string filename, double xsize, double ysize
% Outputs: complex data, raw dimension [no_expts,no_echoes,no_slices,no_views,no_views_2,no_samples], MRD/PPR parameters
% Author: Ruslan Garipov
% Date: 01/03/2014 - swapped views and views2 dimension - now correct
% Date: 01/12/2016 - Ben Hoff (BAH) - streamlined data & par organization

function [img,dim,par] = Get_mrd_3D1_BH(filename)
% Read in wrapped imaged and sensitivity maps
% reordering1, 2 is 'seq' or 'cen'
% reordering1 is for 2D (views)
% reordering2 is for 3D (views2)

% if nargin<2
%     reordering1 = 'seq';
% end
% if nargin<3
%     reordering2 = 'seq';
% end

fid = fopen(filename,'r');
dim = fread(fid,4,'int32');
fseek(fid,18,'bof');
datatype = fread(fid,1, 'uint16');
datatype = dec2hex(datatype);
% fseek(fid,48,'bof');
% scaling = fread(fid,1, 'float32');
% bitsperpixel = fread(fid,1, 'uchar');
fseek(fid,152,'bof');
dimarr = fread(fid,2, 'int32');
% dim5 = val(1);
% dim6 = val(2);
fseek(fid,512,'bof');


if size(datatype,2)>1
    onlydatatype = datatype(2);
    iscomplex = 2;
else
    onlydatatype = datatype(1);
    iscomplex = 1;
end
switch onlydatatype
    case '0' 
        dataformat = 'uchar';
    case '1' 
        dataformat = 'schar';
    case '2' 
        dataformat = 'short';
    case '3' 
        dataformat = 'int16';
    case '4' 
        dataformat = 'int32';
    case '5' 
        dataformat = 'float32';
    case '6' 
        dataformat = 'double';
    otherwise
        dataformat = 'int32';
end

% Read all data at once:
num2read = prod(dim)*prod(dimarr)*iscomplex;
[img, count] = fread(fid,num2read,dataformat);
if (count~=num2read)
    error('We have a problem... count ~= num2read');
end
% sample_filename = char(fread(fid,120,'uchar')');
fseek(fid,120,0);
ppr_text = char(fread(fid,Inf,'uchar')');
fclose(fid);

if iscomplex == 2
    img = complex(img(1:2:count),img(2:2:count));
end

% Reshape data into image-based matrix:
img = permute(reshape(img,[dim',dimarr']),[1,3,2,4,5]);

% % Re-order matrix dimentions for centric acquisitions:
% if strcmp(reordering1,'cen')
%     ind = [ (dim(2)/2+1):dim(3) ; dim(2)/2:-1:1];
%     img = img(:,ind,:,:,:,:);
% end
% if strcmp(reordering2,'cen')
%     ind = [ (dim(3)/2+1):dim(3) ; dim(3)/2:-1:1];
%     img = img(:,:,ind,:,:,:);
% end

% Now parse the PPR parameters:
if numel(ppr_text)>0
    ppr_text = textscan(ppr_text,'%s','delimiter',char(13));

    par = struct('filename',filename);
    for i = 1:length(ppr_text{1})
        str = ppr_text{1}{i};
        if strcmp(str(1),':')
            [field_,rem] = strtok(str(2:end));
            tval = strsplit(rem(2:end),', ');
            switch field_
                case {'VAR_ARRAY'}
                    % Variable, values on multiple lines
                    field_ = tval{1};
                    valen = str2double(tval{2});
                    val = nan(1,valen);
                    tval = cellfun(@str2double,tval(3:end));
                    ct = length(tval);
                    val(1:ct) = tval;
                    ii = 1;
                    while ct<valen
                        tval = strsplit(ppr_text{1}{i+ii},', ');
                        tval(cellfun(@isempty,tval)) = [];
                        tval = cellfun(@str2double,tval);
                        n = length(tval);
                        val(ct+(1:n)) = tval;
                        ct = ct+n;
                        ii = ii+1;
                    end
                case {'OBSERVE_FREQUENCY'}
                    % Fix to single string (ignore any more)
                    val = tval{1};
                case {'SAMPLE_PERIOD','SAMPLE_PERIOD_2'}
                    % Variable, fix to single value
                    field_ = tval{1};
                    val = str2double(tval{2});
                case {'SLICE_THICKNESS','SLICE_SEPARATION'}
                    % Keep only second value
                    val = str2double(tval{3});
                otherwise
                    % Variable, multi-value
                    num = str2double(tval{1});
                    if isnan(num)
                        if length(tval)==1
                            val = tval{1};
                        else
                            % First value is the field name:
                            field_ = tval{1};
                            val = tval(2:end);
                        end
                    else
                        % Start with values:
                        val = tval;
                    end
                    % If all numbers, convert to vector
                    if iscell(val)
                        nval = cellfun(@str2double,val);
                        if ~any(isnan(nval))
                            val = nval;
                        elseif length(val)==1
                            val = val{1};
                        end
                    end
            end
            par.(matlab.lang.makeValidName(field_)) = val;
        end
    end
    if isfield(par,'OBSERVE_FREQUENCY')
        C = textscan(par.OBSERVE_FREQUENCY, '%q');
        par.Nucleus = C{1}{1};
    else
        par.Nucleus = 'Unspecified';
    end
    par.datatype = datatype;
    file_pars = dir(filename);
    par.date = file_pars.date;
else par = [];
end
