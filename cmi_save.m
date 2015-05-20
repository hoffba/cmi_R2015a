function status = cmi_save(maskflag,img,label,fov,fname)
%% Save image data from umich2 GUI as various data types
% Includes FLD, DICOM, TIFF, MASK, MAT

if nargin > 1
    status = false;
    if nargin < 5
        fname = '';
    end
    [pname,name,ext] = fileparts(fname);
    if isempty(pname) % Ask user how to save data
        if maskflag % save mask
            [fname,pname] = uiputfile({'*.mhd','MHD';'*.mask','MASK';'*.fld','FLD';'*.hdr','ANALYZE'},...
                                      'Save MASK as:',name);
            img = img(:,:,:,1);
        else % save image
            [fname,pname] = uiputfile({'*.mhd','MHD';'*.fld','FLD';'*.dcm','DICOM';...
                                       '*.hdr','ANALYZE';'*.vff','VFF';'*.tif','TIF';...
                                       '*.mat','MAT'},...
                                       'Save IMAGE as:',name);
        end
    else
        fname = [name ext];
    end
    if fname
        [~,~,ext] = fileparts(fname);
        if isempty(ext) % Default to save FLD
            fname = [fname '.fld'];
            ext = 'fld';
        else
            ext = ext(2:end);
        end
        % Call format save function
        switch ext
            case 'fld'
                if (nargin < 3)
                    if maskflag
                        label = 'label= "VOI"';
                    else
                        label = 'label= "Image"';
                    end
                end
                if (nargin < 4)
                    fov = size(img);
                    if (length(fov) > 3)
                        fov(4:end) = [];
                    end
                end
                status = saveFLD(fullfile(pname,fname),img,label,fov);
            case 'mask'
                status = saveMASK(fullfile(pname,fname),img);
            case 'vff'
                status = saveVFF(fullfile(pname,fname),img,label,fov);
            case 'mat'
                status = saveMAT(fullfile(pname,fname),img,label,fov);
            case 'tif'
                status = saveTIFF(fullfile(pname,fname),img,label,fov);
            case 'dcm'
                status = saveDICOM(fullfile(pname,fname),img,label,fov);
            case 'mhd'
                status = saveMHD(fullfile(pname,fname),img,label,fov);
            case 'hdr'
                status = saveANALYZE(fullfile(pname,fname),img,label,fov);
        end
        if status
            status = fullfile(pname,fname);
        end
    end

end
