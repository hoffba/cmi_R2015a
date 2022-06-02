function res = job_tPRM(fn_prm,fn_seg,sv_path)
% Perform tPRM processing 
% Inputs: fn_prm = full filename of PRM classification map (as seen by GL)
%         fn_seg = full filename of segmentation map (as seen by GL)
%         sv_path = path to save resulting tPRM maps into
    
    res = [];

    res.ID = extractBefore(flip(extractBefore(flip(fn_seg),filesep)),'.segmentation.nii');

    % Load PRM map (4-color from ImBio SPIROMICS)
    fprintf('Loading PRM map: %s\n',fn_prm);
    info_prm = niftiinfo(fn_prm);
    prm = niftiread(info_prm);
    voxsz = info_prm.PixelDimensions;

    % Load lobe segmentation
    fprintf('Loading segmentation: %s\n',fn_seg);
    info_seg = niftiinfo(fn_seg);
    seg = niftiread(info_seg);

    % Set up labels for results
    info_sv = info_prm;
    info_sv.Datatype = 'int16';
    info_sv.BitsPerPixel = 16;
    prmlabel = [ "Norm" , "fSAD", "emph", "NS" ];
    nprm = numel(prmlabel);
    MF_label = [ "V" , "S" , "B" , "X" ];
    scale=[10000 10000 1e5 1e5];
    
    % Values for 4-color PRM
    T = lobeLoop(seg,@(mask,prm,vals,tag)tabulatePRM(mask,prm,vals,tag),prm,1:4,prmlabel);
    res = lobeTable2struct(T,res,2:5);
        
    % Check if tPRM files already exist
    fnames = string(res.ID) + '.tprm.' + lower(prmlabel)' + '.' + lower(MF_label) + '.nii';
    exist_chk = cellfun(@(x)exist(x,'file'),fullfile(sv_path,fnames+'.gz'));
    
    % Generate MF values
    if ~all(exist_chk,'all')
        fprintf('Calculating Minkowski functionals ...\n');
        p = minkowskiFun(prm,'thresh',1:nprm,'tmode','==','n',[10,10,10],'gridsp',[5,5,5],'voxsz',voxsz,'mask',seg>0);
    end
    
    % Re-grid, scale, and save as .nii.gz
    for iprm = 1:nprm
        
        % tPRM
        for imf = 1:4
            
            fn_out = fnames(iprm,imf);
            svname = fullfile(sv_path,fn_out);
            str = prmlabel(iprm) + '_' + MF_label(imf);
            
            if exist_chk(iprm,imf)
                fprintf('Reading tPRM map from file: %s\n',fn_out);
                try
                    prm = niftiread(svname+'.gz');
                catch
                    fprintf(' ... FAILED.\n');
                    exist_chk(iprm,imf) = false;
                end
            end
            if ~exist_chk(iprm,imf)
                fprintf('Generating tPRM map from grid: %s\n',fn_out);
                prm = grid2img(p.MF(iprm,imf,:),p.ind,seg>0,3,1) * scale(imf);
            end
            
            % Calculate segmentation means:
            T = lobeLoop(seg,@(mask,tprm,str)tabulateTPRM(mask,tprm,str),prm,str);
            res = lobeTable2struct(T,res,2);
            
            % write to /tmpssd or /tmp folder to save read/write time on the zip
            if ~exist_chk(iprm,imf)
                fprintf('Saving MF result: %s.gz\n',svname);
                niftiwrite(int16(prm),svname,info_sv,'Compressed',true);
            end
            
        end
    end
    
