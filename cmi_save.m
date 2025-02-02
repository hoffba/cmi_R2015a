function status = cmi_save(maskflag,img,label,fov,orient,fname)
%% Save image data from umich2 GUI as various data types
% Includes FLD, DICOM, TIFF, MASK, MAT

if nargin > 1
    status = false;
    if nargin < 6
        fname = '';
    end
    
    if maskflag
        opts = struct('prompt','Save MASK as:',...
                      'formats',{{ '*.nii.gz;*.nii', 'NIFTI';...
                                   '*.mhd',    'MHD';...
                                   '*.mask',   'MASK';...
                                   '*.fld',    'FLD';...
                                   '*.hdr',    'ANALYZE'  }});
        img = img(:,:,:,1);
    else
        opts = struct('prompt','Save IMAGE as:',...
                      'formats',{{ '*.nii.gz;*.nii', 'NIFTI';...
                                   '*.mhd',    'MHD';...
                                   '*.fld',    'FLD';...
                                   '*.dcm',    'DICOM';...
                                   '*.hdr',    'ANALYZE';...
                                   '*.vff',    'VFF';...
                                   '*.tif',    'TIF';...
                                   '*.mat',    'MAT'  }});
    end
    
    [fpath,fname,ext] = fileparts(fname);
    if strcmp(ext,'.gz')
        [~,fname,ext2] = fileparts(fname);
        ext = [ext2,ext];
    end
    if isempty(fpath) % Ask user how to save data
        [fname,fpath,indx] = uiputfile(opts.formats,opts.prompt,fname);
    else
        fname = [fname,ext];
        indx = find(cellfun(@(x)contains(x,ext),opts.formats(:,1)),1);
    end
    if fname
        if indx
            fmt = opts.formats{indx,2};
        else % Default to save NIfTI
            fname = [fname '.nii.gz'];
            fmt = 'NIFTI';
        end
        fname = fullfile(fpath,fname);
        
        if nargin<3
            if maskflag
                label = 'VOI';
            else
                label = 'Image';
            end
            if strcmp(fmt,'FLD')
                label = ['label= "',label,'"'];
            end
        end
        if nargin<4
            fov = size(img);
            if (length(fov) > 3)
                fov(4:end) = [];
            end
        end
        
        % Call format save function
        status = feval(['save',fmt],fname,img,label,fov,orient);
        if status
            status = fullfile(fpath,fname);
        end
    end

end
