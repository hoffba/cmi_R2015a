function Mreg = pipeline_inverseTform(procdir,fn_img,info_hom,flag_nn)
% Transform images from Ref space (Exp) to Homologous space (Ins)
% Inputs:
%   procdir  = case processing directory
%   fn_img   = image file to transform (Exp space)
%   info_hom = image Nifti info (Ins space)
%   flag_nn  = T/F flag for nearest neighbor

Mreg = [];
if nargin<5
    flag_nn = true;
end

[~,ID] = fileparts(procdir);
fn_ref = fullfile(procdir,[ID,'.exp.nii.gz']);
fn_seg = fullfile(procdir,[ID,'.exp.label.nii.gz']);
path_elx = fullfile(procdir,['elastix_',ID]);

if isfile(fn_ref) && isfile(fn_seg) && isfolder(path_elx)

    % Inverse transform
    fn_tfi = dir(fullfile(path_elx,'InverseTransformParameters.*.txt'));
    if isempty(fn_tfi)
        fn_tf = dir(fullfile(path_elx,'TransformParameters.*.txt'));
        if isempty(fn_tf)
            fprintf('No transform found in: %s\n',path_elx);
            fn_tfi = '';
        else
            fprintf('Generating inverse transform.\n');
            fn_tfi = inverseTransform(fullfile(path_elx,fn_tf(end).name),fn_ref,fn_seg,info_hom);
        end
    else
        fprintf('InverseTransformParameters.*.txt files found.\n');
        fn_tfi = fullfile(path_elx,fn_tfi(end).name);
    end
    
    if ~isempty(fn_tfi)
        fprintf('Transforming Ref->Hom ... ');
        Mreg = transformImage(fn_tfi,flag_nn,fn_img);
        fprintf('done\n');
    end

end