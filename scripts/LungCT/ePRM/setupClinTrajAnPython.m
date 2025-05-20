function python_path = setupClinTrajAnPython()
% Sets up a Python environment for use with ClinTrajAn

% Path to ClinTrajAn Python environment
cta_env_path = 'C:\Users\bahoff\.conda\envs\ClinTrajAn';

if isfolder(cta_env_path)
    % Return existing environment path
    python_path = fullfile(cta_env_path,'python.exe');
else
    % Need to create new environment for ClinTrajAn
    
end
