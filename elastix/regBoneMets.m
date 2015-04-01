% BH script
function stat = regBoneMets(cmiObj0,cmiObj1,waitchk)
%   Input: cmiObj0 = CMIclass object containing Expiration image
%          cmiObj1 = CMIclass object containing Inspiration image

stat = false;
if (nargin>=2) && isa(cmiObj0,'CMIclass') && isa(cmiObj1,'CMIclass') ...
        && cmiObj0.img.check && cmiObj0.img.mask.check ...
        && cmiObj1.img.check && cmiObj1.img.mask.check
    
    waitstr = '''';
    if (nargin<3) || ~waitchk
        waitstr = ';csh''&';
    end
    
    % Determine filenames
    outdir = cmiObj0.img.dir;
    [~,dstr] = fileparts(outdir);
    if (length(dstr)>2) && any(strcmpi(dstr(end-2:end),{'exp','ins'}))
        outdir = fileparts(outdir);
    end
    outdir = fullfile(outdir,[cmiObj1.img.name,'_elxreg']);
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end
    
    % Temporary files saved as MHD for Elastix compatibility
    ext = '.mhd';
    homfn = fullfile(outdir,['elxtemp-mImg',ext]);
    reffn = fullfile(outdir,['elxtemp-fImg',ext]);
    fmaskfn = fullfile(outdir,['elxtemp-fMask',ext]);
    mmaskfn = fullfile(outdir,['elxtemp-mMask',ext]);
    origmmaskfn = fullfile(outdir,['elxtemp-origmMask',ext]);

    % Temporarily save .mhd images
    fov0 = cmiObj0.img.voxsz.*cmiObj0.img.dims(1:3);
    fov1 = cmiObj1.img.voxsz.*cmiObj1.img.dims(1:3);
    disp('Saving temporary image files ...');
    % Un-clamped Homologous image
    str = cmi_save(0,cmiObj1.img.mat(:,:,:,cmiObj1.vec),{'Hom'},fov1,homfn);
    [~,str] = fileparts(str);
    disp(['Saved ',str]);
    % Un-clamped Reference image
    str = cmi_save(0,cmiObj0.img.mat(:,:,:,cmiObj0.vec),{'Ref'},fov0,reffn);
    [~,str] = fileparts(str);
    disp(['Saved ',str]);
    % VOI
    str = cmi_save(1,cmiObj1.img.mask.mat(:,:,:,1)*1000,{'VOI'},fov1,origmmaskfn);
    [~,str] = fileparts(str);
    disp(['Saved ',str]);
    cmiObj0.img.mask.morph('dilate',[5 2*(cmiObj0.img.dims(3)>20)]);
    str = cmi_save(1,cmiObj0.img.mask.mat(:,:,:,1)*1000,{'VOI'},fov0,fmaskfn);
    [~,str] = fileparts(str);
    disp(['Saved ',str]);
    cmiObj1.img.mask.morph('dilate',[5 2*(cmiObj1.img.dims(3)>20)]);
    str = cmi_save(1,cmiObj1.img.mask.mat(:,:,:,1)*1000,{'VOI'},fov1,mmaskfn);
    [~,str] = fileparts(str);
    disp(['Saved ',str]);

    % Determine initial translation transform based on VOI center of mass:
    [xx,yy,zz] = ind2sub(cmiObj0.img.dims(1:3),find(cmiObj0.img.mask.mat));
    com0 = ([mean(xx),mean(yy),mean(zz)]-cmiObj0.img.dims(1:3)/2).*cmiObj0.img.voxsz;
    [xx,yy,zz] = ind2sub(cmiObj1.img.dims(1:3),find(cmiObj1.img.mask.mat));
    com1 = ([mean(xx),mean(yy),mean(zz)]-cmiObj1.img.dims(1:3)/2).*cmiObj1.img.voxsz;
    t = com1 - com0;
    initT = fullfile(outdir,'InitialTransform.txt');
    saveAffineTxformPars(t,struct('dims',cmiObj0.img.dims(1:3),...
                                  'voxsz',cmiObj0.img.voxsz),...
                         initT);
    
    % Create Elastix Parameter File:
    p = initAffineParStruct; % Initial settings are for mice
%     if max(cmiObj0.img.voxsz/min(cmiObj0.img.voxsz))>2
        p.Pyramid.FixedImagePyramidSchedule = ones(1,p.Pyramid.NumberOfResolutions*3);
%     end
%     if max(cmiObj1.img.voxsz/min(cmiObj0.img.voxsz))>2
        p.Pyramid.MovingImagePyramidSchedule = ones(1,p.Pyramid.NumberOfResolutions*3);
%     end
%     p.Outputs.WriteResultImageAfterEachResolution = 'true';
    pstr = saveElxParFile(p,fullfile(outdir,'ElastixParameters.txt'));
%     pstr = '/mnt/cmi/projects/CMI/CMI_Matlab_Programs/BenHoff/cmi/elastix/BoneMets/ClinBoneMets.affine.txt';

    % Window Name:
    namestr = ['Elastix Registration: ',cmiObj1.img.name,' --> ',cmiObj0.img.name];
    % Call to Elastix:
    elxstr = ['/opt/elastix/bin/elastix',...
              ' -f ',reffn,...
              ' -m ',homfn,...
              ' -fMask ',fmaskfn,...
              ' -mMask ',mmaskfn,...
              ' -out ',outdir,...
              ' -t0 ',initT,...
              ' -p ',pstr,'; '];
    
    % Executes in new xterm window
    %  -geometry sets window size
    %  -T sets window name
    %  -e executes commands (within ' '), must be last xterm input
    %  csh at end to keep window open
    %  & at end to free Matlab from waiting
    stat = system(['xterm -geometry 170x50 -T "',namestr,'"',...
                   ' -e ''',elxstr,waitstr]);
               
    % VOI transform needs nearest neighbor interpolation
    fid = fopen(fullfile(outdir,'TransformParameters.0.txt'),'r');
    str = fread(fid,'*char')';
    fclose(fid);
    ind0 = strfind(str,'(ResampleInterpolator ')+22;
    ind1 = strfind(str(ind0+1:end),'"')+ind0;
    str = [str(1:ind0),'FinalNearestNeighborInterpolator',str(ind1:end)];
    fid = fopen(fullfile(outdir,'TransformParametersNN.0.txt'),'w');
    fwrite(fid,str,'char');
    fclose(fid);
    % Call to Transformix: Transform moving VOI
    tfxstr = ['/opt/elastix/bin/transformix',...
              ' -out ',outdir,...
              ' -in ',origmmaskfn,...
              ' -tp ',fullfile(outdir,'TransformParametersNN.0.txt'),'; '];
    % Cleanup string (removes temporary image files):
    custr = ['find ',outdir,' -name "elxtemp-*" -exec rm -f {} \; '];
    stat = system(['xterm -geometry 170x50 -T "',namestr,'"',...
                   ' -e ''',tfxstr,custr,waitstr]);
end
