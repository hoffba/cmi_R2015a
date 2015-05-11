function status = saveTIFF(fname,img,label,~,~)

status = true;
[fpath,bname,~] = fileparts(fname);

mmax = squeeze(max(max(max(img))));
[d(1),d(2),d(3),d(4)] = size(img);

for j = 1:d(4)
    
    % TIFF cannot be 4D, so separate and save with label in name:
    oname = fullfile(fpath,[bname,'_',label{j},'.tif']);
    
    % Determine needed bitdepth:
    bd = 'uint16';
    chk = find([ (mmax(j)<1) , (mmax(j)>(2^16)) ],1);
    if ~isempty(chk)
        answer = questdlg('Values are not scaled well for TIFF, auto-scale?',...
            'TIFF scaling','Yes','No','Cancel','Yes');
        switch answer
            case 'Cancel'
                break;
            case 'Yes'
                img(:,:,:,j) = img(:,:,:,j) * ((2^16)-1) / mmax(j);
            case 'No'
                if chk==1
                    bd = 'uint8';
                end
        end
    elseif mmax<256
        bd = 'uint8';
    end
    
    for i = 1:d(3)
        imwrite(cast(img(:,:,i,j),bd),oname,...
                'Compression','none',...
                'WriteMode','append');
    end
end