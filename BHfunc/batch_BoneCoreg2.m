function batch_BoneCoreg2(cmiRef,cmiHom)
% All images and VOIs should already be saved as MHD format

% Select all images to coregister to baseline
[fnames,bdir] = uigetfile('*.mhd','Select images to coregister:',...
                            '','MultiSelect','on');
                        
if iscellstr(fnames)
    
    % Load fixed image and mask:
    fixedfn = fullfile(bdir,fnames{1});
    [~,ffnbase] = fileparts(fixedfn);
    cmiRef.loadImg(false,fixedfn);
    cmiRef.loadMask(fullfile(bdir,[ffnbase,'_VOI.mhd']));
    
    for ifn = 2:length(fnames)
        
        % Load moving image and mask:
        [~,mfnbase] = fileparts(fnames{ifn});
        fprintf('\nProcessing %s\n',mfnbase);
        cmiHom.loadImg(false,fullfile(bdir,fnames{ifn}));
        cmiHom.loadMask(fullfile(bdir,[mfnbase,'_VOI.mhd']));
        
        % Perform coregistration:
        % Files saved into [ fname '_elxreg' ]
        elxdir = fullfile(bdir,[cmiHom.img.name,'_elxreg']);
        regBoneMets2(cmiRef,cmiHom);
        
        if exist(fullfile(elxdir,'result.0.mhd'),'file')
            % Save coregistered image and VOI back in base directory:
            result_img = fullfile(elxdir,'result.0.mhd');
            result_voi = fullfile(elxdir,'result.mhd');
            cmiRef.loadImg(false,result_img);
            cmiRef.loadMask(result_voi);
            fixedfn = fullfile(bdir,[mfnbase,'_RR.mhd']);
            cmiRef.saveImg(fixedfn);
            cmiRef.img.saveMask(fullfile(bdir,[mfnbase,'_RR_VOI.mhd']));
            delete(result_img,result_voi,[result_img(1:end-3),'raw'],[result_voi(1:end-3),'raw']);
        else
            disp('Coregistration failed ... stopping batch process');
            break;
        end
    end
end