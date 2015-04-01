function status = saveVFF(fname,img,label,max_ext)
%% Save image as VFF format, compatible with GEHC MicroView

d = size(img);
if length(d)>3
    nv = d(4);
    d = d(1:3);
else
    nv = 1;
end

if length(d)==3 && length(max_ext)==3 && ~isempty(fname)
    dres = 16; % format bit depth
    voxsz = max_ext ./ d;
    [pathstr, name, ext] = fileparts(fname);
    % add VFF file extension to fname if necessary
    if (isempty(ext) || ~strcmp(ext,'.vff'))
        fname = [pathstr name '.vff'];
    end
    fido = fopen(fname,'w','b');
    % Write VFF header information
    fprintf(fido,'%s\n',...
        'ncaa',...
        'rank=3;',...
        'type=raster;',...
        'format=slice;',...
        sprintf('bands=%i;',nv),...
        sprintf('bits=%i',dres),...
        sprintf('size=%i %i %i;',d),...
        'origin=0 0 0;',... % not sure how to put this yet
        sprintf('spacing=%f %f %f;',voxsz),...
        'elementsize=1.00;',...
        'modality=4D;',...
        ['title=' sprintf(' "%s"',label{:}) ';'],...
        'description=Saved from CMI Matlab code (University of Michigan);',...
        char(12)); % end of header marker
    % Write VFF image data to file
    status = fseek(fido,0,'eof');
    if (status==-1)
        ferror(fido);
    else
        if nv>1
            img = shiftdim(img,1);
        end
        fwrite(fido,img,'short');
    end
    fclose(fido);
end
