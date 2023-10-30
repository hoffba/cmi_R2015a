function seg = brainSeg_SPM(procdir,fn_img,info)
seg = [];

% First check that you have SPM loaded and mapped
spmloc = which('spm');
if isempty(spmloc)
    warning('SPM is not available on this system.');
    return;
end

if ischar(fn_img)
    if isfolder(fn_img)
        % DICOM

    elseif isfile(fn_img)
        
    end
elseif isnumeric(fn_img) && nargin==2
else
    fprintf('\nInvalid inputs.\n')
    return;
end

% Run the job in SPM
jobfile = {writeJobFile(spmloc,fn_img,procdir)};
spm('defaults','FMRI');
spm_jobman('run', jobfile);


function fn = writeJobFile(spmloc,fn_img,procdir)
fn_img = regexprep(fn_img,'\\','\\\\');
fn = fullfile(procdir,'SPMjob.m');
fid = fopen(fn,'w');
fprintf(fid,['matlabbatch{1}.spm.spatial.preproc.channel.vols = {''%s,1''};\n',...
             'matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;\n',...
             'matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;\n',...
             'matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];\n'],fn_img);
for i = 1:6
    fntpm = fullfile(fileparts(spmloc),'tpm','TPM.nii');
    fntpm = regexprep(fntpm,'\\','\\\\');
    fprintf(fid,['matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {''%s,%u''};\n',...
                 'matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;\n',...
                 'matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];\n',...
                 'matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];\n'],fntpm,i);
end
fprintf(fid,['matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;\n',...
             'matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;\n',...
             'matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];\n',...
             'matlabbatch{1}.spm.spatial.preproc.warp.affreg = ''mni'';\n',...
             'matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;\n',...
             'matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;\n',...
             'matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];\n']);
fclose(fid);