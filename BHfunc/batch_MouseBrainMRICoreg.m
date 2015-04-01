function batch_MouseBrainMRICoreg(cmiRef,cmiHom)
% All images and VOIs should already be saved as MHD format

% Select all images to coregister to baseline
[fnames,bdir] = uigetfile('*.mhd','Select images to coregister:',...
                            '','MultiSelect','on');
                        
if iscellstr(fnames)
    
    % Load fixed image and mask:
    fixedfn = fullfile(bdir,fnames{1});
    [~,ffnbase] = fileparts(fixedfn);
    fixedmaskfn = fullfile(bdir,[ffnbase,'_VOI.mhd']);
    cmiRef.loadImg(false,fixedfn);
    cmiRef.loadMask(fixedmaskfn);
    cmiRef.img.mask.morph('dilate',[5,0]);
    cmiRef.img.interp([256,256,cmiRef.img.dims(3)]);
    
    for ifn = 2:length(fnames)
        
        % Load moving image and mask:
        [~,mfnbase] = fileparts(fnames{ifn});
        fprintf('\nProcessing %s\n',mfnbase);
        cmiHom.loadImg(false,fullfile(bdir,fnames{ifn}));
        cmiHom.img.interp([256,256,cmiHom.img.dims(3)]);
        
        % Perform coregistration:
        % Files saved into [ fname '_elxreg' ]
        elxdir = fullfile(bdir,[cmiHom.img.name,'_elxreg']);
        regMouseBrainMRI(cmiRef,cmiHom,true);
        
        if exist(fullfile(elxdir,'result.mhd'),'file')
            % Save coregistered image and VOI back in base directory:
            result_img = fullfile(elxdir,'result.mhd');
            cmiRef.loadImg(false,result_img);
            cmiRef.loadMask(fixedmaskfn);
            cmiRef.img.mask.morph('dilate',[5,0]);
            fixedfn = fullfile(bdir,[mfnbase,'_RR.mhd']);
            cmiRef.saveImg(fixedfn);
            delete(result_img,[result_img(1:end-3),'raw']);
        else
            disp('Coregistration failed ... stopping batch process');
            break;
        end
    end
end