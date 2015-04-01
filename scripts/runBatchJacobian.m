% CMI script
function runBatchJacobian(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings
% 
% Revision history:
% 10 july 2014 jlb created to read excel file of names and create script 
%                  file to run MiamiFUSE to generate Jacobian

%% Load  names plus additional data from an Excel file
[fnames, exp_mean_blood, exp_mean_air, ins_mean_blood, ins_mean_air] = read_excel_study_to_PRM_fullfnames();
assignin('base','fnames',fnames);  % for debugging only

if isempty(fnames)
    fprintf('ERROR: no filenames to process from excel file. Exiting script now.\n');
    return;
end

%% Create script file

% cycle through all fnames
for i=1:size(fnames,2)
    % Open file
    out_tname = 'script_tmp.cli';
    out_tpath = pwd;
    sfname= fullfile(out_tpath,out_tname);
    fid = fopen(sfname,'w');
    if fid==-1
        fprintf('****ERROR opening %s\n',out_tname);
        return;
    end
    
    % Create script file
    fprintf(fid,'net_clear\n');
    fprintf(fid,'net_read /mnt/cmi/projects/CT_Lung/clinical/COPDGene/network/Jacobian_simple.net\n');
    fprintf(fid,'AVSmessage Ok\n');
    fprintf(fid,'parm_set "read field.user.11":"Read Field Browser" %s\n',fnames{i}{2});
    fprintf(fid,'parm_set "point cluster(s).user.2":"Read dat Feature File" %s\n',fnames{i}{3});
    fprintf(fid,'parm_set "point cluster(s).user.2":"Read dat Feature File" %s\n',fnames{i}{4});
    fprintf(fid,'AVSmessage Ok\n');
    fprintf(fid,'parm_set "point cluster(s).user.1":"Read dat Feature File" %s\n',fnames{i}{5});
    fprintf(fid,'parm_set "point cluster(s).user.1":"Read dat Feature File" %s\n',fnames{i}{6});
    fprintf(fid,'AVSmessage Ok\n');
    fprintf(fid,'parm_set oneshot.user.7:oneshot 1\n');
    out_pname = fnames{i}{1};
    out_ppath = out_pname;
    out_pname = out_pname(end-5:end);
    out_pname = strcat(out_pname,'_Exp_Jac.fld');
    fprintf(fid,'parm_set "write field.user.16":"Write Field Browser" %s\n',fullfile(out_ppath,out_pname));
    %fprintf(fid,'net_clear\n');
    fclose(fid);
    
    % Run script file
    unix_command = sprintf('miamiscript -cli "script -play %s -quit"',sfname);
    status = unix(unix_command);

end