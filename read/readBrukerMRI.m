function [img,label,fov,p] = readBrukerMRI(varargin)

imname = varargin{1};
[fdir,tname] = fileparts(imname);

if strcmp(tname,'2dseq')
    % Recontructed image
    flag = true;
    p = readpar(fullfile(fileparts(fileparts(fdir)),'method'));
    pr = readpar(fullfile(fdir,'reco'));
    p.d_reco = pr.d;
    d = pr.d;
    fov = pr.fov*10;
    if length(d)==2
        d(3) = p.ns;
        fov(3) = p.thk;
    end
    d(4) = pr.nrep * pr.ns;
elseif exist(fullfile(fdir,'method'),'file')
    % Raw k-space data
    flag = false;
    p = readpar(fullfile(fdir,'method'));
    d = p.d;
    fov = p.fov;
    if length(d)==2
        d(3) = p.ns;
        fov(3) = p.thk*p.ns;
    end
    d(4) = p.nrep;
else
    error('Invalid directory input.');
end
    
% Read binary data file:
fid = fopen(imname,'r');
img = fread(fid,'int16');
fclose(fid);

% Adjust image based on method/reco parameters:
label = strcat(p.seq(9:end-1),cellfun(@num2str,num2cell(1:d(4)),'UniformOutput',false));
img = reshape(img,[d(2),d(1),d(3),d(4)]);
if flag
    for i = 1:d(4)
        img(:,:,:,i) = img(:,:,:,i) / pr.scale(i) - pr.offset(i);
    end
end


% Subfunction for reading parameters from file:
function p = readpar(fname)
if exist(fname,'file')
    [~,fstr] = fileparts(fname);
    switch fstr
        case 'method'
            pnames = {'PVM_Matrix'          ,'d';...
                      'PVM_NRepetitions'    ,'nrep';...
                      'PVM_EffSWh'          ,'sw';... % Hz
                      'PVM_NEchoImages'     ,'ne';...
                      'PVM_Fov'             ,'fov';... % mm
                      'PVM_SliceThick'      ,'thk';... % mm
                      'PVM_EchoTime'        ,'te';...
                      'EchoSpacing'         ,'te2';...
                      'PVM_RepetitionTime'  ,'tr';...
                      'PVM_NSPacks'         ,'nsp';...
                      'PVM_SPackArrNSlices' ,'ns';...
                      'PVM_SPackArrSliceGap','gap';...
                      'PVM_DummyScans'      ,'ndum';...
                      'PVM_FatSupOnOff'     ,'fatsup';...
                      'PVM_MagTransOnOff'   ,'mt';...
                      'PVM_MagTransOffset'  ,'mt_off';...
                      'PVM_MagTransPower'   ,'mt_pow';...
                      'PVM_MagTransPulsNumb','mt_n';...
                      'PVM_MagTransInterDelay','mt_delay';
                      'Method'              ,'seq'};
        case 'reco'
            pnames = {'RECO_size'               ,'d';...
                      'RECO_sw'                 ,'sw';... % Hz
                      'RecoNumRepetitions'      ,'nrep';...
                      'RecoObjectsPerSetupStep' ,'ns';...
                      'RECO_fov'                ,'fov';...
                      'RECO_map_slope'          ,'scale';...
                      'RECO_map_offset'         ,'offset'};
        case 'acqp'
    end
    go = true;
    pflag = false(size(pnames,1),1);
    fid = fopen(fname,'r');
    while go
        line = fgetl(fid);
        val = regexp(line,'##\$(.*)=(.*)','tokens');
        if ~isempty(val)
            ind = find(strcmp(val{1}{1},pnames(:,1)),1);
            if ~isempty(ind)
                val = val{1}{2};
                sval = str2double(val);
                if ~isnan(sval)
                    val = sval;
                elseif strcmp(val(1:2),'( ')
                    n = sscanf(val(3:end-2),'%f, %f, %f, %f')';
                    if length(n)==1
                        n = [1,n];
                    end
                    val = nan(n);
                    ct = 0;
                    while ct<prod(n)
                        tval = sscanf(fgetl(fid),'%f');
                        nval = length(tval);
                        val(ct+(1:nval)) = tval;
                        ct = ct+nval;
                    end
                end
                p.(pnames{ind,2}) = val;
                pflag(ind) = true;
            end
        end
        if feof(fid) || all(pflag)
            go = false;
        end
    end
    fclose(fid);
else
    error('File not found: %s',fname)
end


