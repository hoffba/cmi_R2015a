function stat = write_sbatch(fname,jobname,username,cores,mem,fcn,inputs_fname)

stat = false;

GL_turbo_path = '/nfs/turbo/umms-cgalban';
GL_Matlab_path = [GL_turbo_path,'/GreatLakes/cmi_R2015a'];

% Check that inputs_fname is on Turbo and fix path if necessary
if ~startsWith(inputs_fname,GL_turbo_path)
    errstr = '';
    tokens = regexp(inputs_fname,'([A-Z]+:)(.*)','tokens');
    if ~(isempty(tokens) || isempty(tokens{1}))
        driveLetter = tokens{1}{1};
        fpath = tokens{1}{2};
        netdrives = findNetDrives;
        ind = strcmp(driveLetter,{netdrives.Drive});
        if contains(netdrives(ind).Remote,'umms-cgalban')
            inputs_fname = [GL_turbo_path,regexprep(fpath,'\\','/')];
        else
            errstr = 'inputs_fname must be located in Turbo storage.';
        end
    else
        errstr = 'Invalid path for inputs_fname';
    end
    if ~isempty(errstr)
        error(errstr);
    end
end

fprintf('Writing SBATCH file: %s\n',inputs_fname);

str = [ '#!/bin/bash\n\n',... % The interpreter used to execute the script
        '#“#SBATCH” directives that convey submission options:\n\n',...
        sprintf('#SBATCH --job-name=%s\n',jobname),...
        '#SBATCH --account=cgalban1\n',...
        '#SBATCH --partition=standard\n',...
        '#SBATCH --nodes=1\n',...
        '#SBATCH --ntasks-per-node=1\n',...
        sprintf('#SBATCH --cpus-per-task=%u\n',min(cores,36)),... % Standard node has 36 cores
        sprintf('#SBATCH --mail-user=%s@med.umich.edu\n',username),...
        '#SBATCH --mail-type=END\n',...
        sprintf('#SBATCH --mem=%ug\n',min(mem,180)),... % Standard node has max 180G requestable
        '#SBATCH --output=./%%x-%%j\n',...
        '\n',...
        'if [[ $SLURM_JOB_NODELIST ]] ; then\n',...
        '   echo "Running on"\n',...
        '   scontrol show hostnames $SLURM_JOB_NODELIST\n',...
        '   cat coinflip.sbat\n',...
        'fi\n\n',...
        'module load matlab\n\n',...
        'my_job_header\n\n',...
        'echo -n "My job array ID is:  "\n',...
        'env | grep ARRAY_TASK_ID\n\n',...
        '#  Put your job commands after this line\n',...
        sprintf(['matlab -nodisplay -r "cd(''%s'');',...
                    'cmi_setPath;%s(''%s''); exit"\n'],...
                    GL_Matlab_path,fcn,inputs_fname)];
fid = fopen(fname,'w');
if fid
    fprintf(fid,str);
    fclose(fid);
    stat = true;
end
