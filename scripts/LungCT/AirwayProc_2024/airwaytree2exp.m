function xyz = airwaytree2exp(procdir)
% Transforms airway nodes (real and simulated) to EXP space
% ** Note: Node coordinates are in non-oriented units
% Inputs:
%   procdir = path of processing case

[~,ID] = fileparts(procdir);
fn_tree = fullfile(procdir,[ID,'.AirwayProc'],[ID,'_AirwayTreeSim.mat']);
if isfile(fn_tree)
    p = load(fn_tree);
    if isfield(p,'N_exp')
        xyz = p.N_exp(:,2:4);
    else
        % Filenames relevant to registration
        fn_ref = fullfile(procdir,[ID,'.exp.nii.gz']);
        fn_seg = fullfile(procdir,[ID,'.exp.label.nii.gz']);
        fn_hom = fullfile(procdir,[ID,'.ins.nii.gz']);
        path_elx = fullfile(procdir,['elastix_',ID]);

        if isfolder(path_elx)
            % Find inverse transform
            fn_tfi = dir(fullfile(path_elx,'InverseTransformParameters.*.txt'));
            if isempty(fn_tfi)
                fn_tf = dir(fullfile(path_elx,'TransformParameteres.*.txt'));
                if isempty(fn_tf)
                    fn_tfi = '';
                else
                    fn_tfi = inverseTransform(fn_tf,fn_ref,fn_seg);
                end
            else
                fn_tfi = fullfile(path_elx,fn_tfi(end).name);
            end
        
            % Transform points with transformix
            if ~isempty(fn_tfi)
                % Find INS node locations to transform
                p = load(fn_tree);
                info = niftiinfo(fn_hom);
                % [~,~,fov,orient,info] = readNIFTI(fn_hom);
                % orient = diag([-1 -1 1 1]) * orient;
                % voxsz = fov ./ info.ImageSize;
                xyz = (orient * [p.N(:,[3,2,4])./info.PixelDimensions,...
                                                ones(size(p.N,1),1)]')';
                N_exp = transformPoints(fn_tfi,xyz(:,1:3));
                N_exp = [p.N(:,1),N_exp(:,[2,1,3])];

                % Save points to airway tree MAT file
                save(fn_tree,'N_exp','-append');
            end
        end

    end
end
