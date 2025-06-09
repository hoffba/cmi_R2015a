function [score, score_case, coeff, explained] = PCApython(savepath, Xs_combined, your_data_size, numComp)
% PCApython - Apply PCA using Python's sklearn implementation with saved parameters

% Check for python setup
if ispref('Clintrajan','PythonPath')
    python_path = getpref('Clintrajan','PythonPath');
else
    python_path = update_python_paths();
end

% Create a temporary directory for data exchange with Python
% Save the data to a temporary MAT file
input_file = fullfile(savepath, 'temp_pca_data.mat');
output_file = fullfile(savepath, 'pca_results.mat');

% Get the ePRM path (same as your pipeline)
path_ePRM = fullfile(fileparts(which('cmi')),'scripts','LungCT','ePRM');

save(input_file, 'Xs_combined', 'your_data_size', 'numComp', 'path_ePRM');

% Write the Python script to a file line by line
script_file = fullfile(savepath, 'run_pca.py');
fid = fopen(script_file, 'w');

% Write each line separately to avoid escape character issues
fprintf(fid, 'import numpy as np\n');
fprintf(fid, 'import pickle\n');
fprintf(fid, 'from sklearn.decomposition import PCA\n');
fprintf(fid, 'import scipy.io as sio\n');
fprintf(fid, 'import sys\n');
fprintf(fid, 'import os\n');
fprintf(fid, '\n');
fprintf(fid, 'def run_pca(input_file, output_file):\n');
fprintf(fid, '    # Load data from MATLAB workspace\n');
fprintf(fid, '    data = sio.loadmat(input_file)\n');
fprintf(fid, '    \n');
fprintf(fid, '    # Extract data\n');
fprintf(fid, '    Xs = data["Xs_combined"]\n');
fprintf(fid, '    your_data_size = int(data["your_data_size"][0,0])\n');
fprintf(fid, '    numComp = int(data["numComp"][0,0])\n');
fprintf(fid, '    pca_model_path = str(data["path_ePRM"][0])\n');
fprintf(fid, '    \n');
fprintf(fid, '    # Print shape for debugging\n');
fprintf(fid, '    print(f"Input data shape: {Xs.shape}, Cases: {your_data_size}, Components: {numComp}")\n');
fprintf(fid, '    print(f"Looking for PCA model in: {pca_model_path}")\n');
fprintf(fid, '    \n');
fprintf(fid, '    # Load saved PCA model (REQUIRED)\n');
fprintf(fid, '    pca_model_file = os.path.join(pca_model_path, "pca_model_complete.pkl")\n');
fprintf(fid, '    pca_params_file = os.path.join(pca_model_path, "pca_parameters.pkl")\n');
fprintf(fid, '    \n');
fprintf(fid, '    # Load files directly - will crash if files don''t exist\n');
fprintf(fid, '    with open(pca_model_file, "rb") as f:\n');
fprintf(fid, '        pca = pickle.load(f)\n');
fprintf(fid, '    \n');
fprintf(fid, '    with open(pca_params_file, "rb") as f:\n');
fprintf(fid, '        pca_params = pickle.load(f)\n');
fprintf(fid, '    \n');
fprintf(fid, '    # Apply saved PCA transformation (no fitting!)\n');
fprintf(fid, '    Y_full = pca.transform(Xs)\n');
fprintf(fid, '    print(f"PCA model has {pca.n_components_} components")\n');
fprintf(fid, '    print(f"Pipeline expects {numComp} components")\n');
fprintf(fid, '    \n');
fprintf(fid, '    # Truncate to the number of components expected by pipeline\n');
fprintf(fid, '    if Y_full.shape[1] > numComp:\n');
fprintf(fid, '        print(f"Truncating from {Y_full.shape[1]} to {numComp} components")\n');
fprintf(fid, '        Y = Y_full[:, :numComp]\n');
fprintf(fid, '        v = pca.components_[:numComp, :].T\n');
fprintf(fid, '        explained = pca.explained_variance_ratio_[:numComp] * 100\n');
fprintf(fid, '    else:\n');
fprintf(fid, '        Y = Y_full\n');
fprintf(fid, '        v = pca.components_.T\n');
fprintf(fid, '        explained = pca.explained_variance_ratio_ * 100\n');
fprintf(fid, '    \n');
fprintf(fid, '    print(f"Using {Y.shape[1]} components for analysis")\n');
fprintf(fid, '    \n');
fprintf(fid, '    print("Using SAVED PCA model - consistent results!")\n');
fprintf(fid, '    \n');
fprintf(fid, '    # Separate case points from reference points\n');
fprintf(fid, '    Y_case = Y[:your_data_size]\n');
fprintf(fid, '    \n');
fprintf(fid, '    print(f"PCA completed successfully with {numComp} components.")\n');
fprintf(fid, '    \n');
fprintf(fid, '    # Save results back to MATLAB\n');
fprintf(fid, '    sio.savemat(output_file, {\n');
fprintf(fid, '        "score": Y,\n');
fprintf(fid, '        "score_case": Y_case,\n');
fprintf(fid, '        "coeff": v,\n');
fprintf(fid, '        "explained": explained\n');
fprintf(fid, '    })\n');
fprintf(fid, '    print(f"Results saved to {output_file}")\n');
fprintf(fid, '\n');
fprintf(fid, 'if __name__ == "__main__":\n');
fprintf(fid, '    if len(sys.argv) != 3:\n');
fprintf(fid, '        print("Usage: python run_pca.py input_file output_file")\n');
fprintf(fid, '        sys.exit(1)\n');
fprintf(fid, '    run_pca(sys.argv[1], sys.argv[2])\n');
fclose(fid);

% Execute Python script with command line arguments
disp('Running PCA using SAVED model...');
command = sprintf('"%s" "%s" "%s" "%s"', ...
    python_path, script_file, input_file, output_file);
[status, commandOutput] = system(command);

% Display output for debugging
disp(commandOutput);
if status ~= 0
    fprintf('Error running Python script:\n%s\n', commandOutput);
    error('Python script failed to execute properly');
end

% Load the results from Python
if exist(output_file, 'file')
    results = load(output_file);
    score = results.score;
    score_case = results.score_case;
    coeff = results.coeff;
    explained = results.explained;
   
    disp('Successfully loaded PCA results from SAVED MODEL');
    fprintf('Final output dimensions: score [%d, %d], score_case [%d, %d]\n', ...
        size(score), size(score_case));
else
    error('PCA results file not found. Python script may have failed.');
end
%Reduction is done
% Clean up temporary files (optional - uncomment if needed)
% delete(input_file);
% delete(output_file);
% delete(script_file);
end