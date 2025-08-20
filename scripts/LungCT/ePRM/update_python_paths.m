function pypath = update_python_paths(x)
    % UPDATE_PYTHON_PATHS Find, store, and retrieve Python path for Clintrajan
    %
    % Usage:
    %   python_path = update_python_paths()     % Get existing path or set new one
    %   python_path = update_python_paths(true) % Force update path
    %   
    % Returns the Python path to use with Clintrajan functions
    
    pypath = '';
    newpath = false;
    force_flag = false;

    % Parse inputs
    if nargin
        if ischar(x)
            if isfile(x)
                pypath = x;
                newpath = true;
            else
                error('Python path does not exist: %s\n',x);
            end
        elseif islogical(x)
            force_flag = x;
        end
    end

    % Check for existing python path
    if isempty(pypath) && ispref('Clintrajan','PythonPath')
        pypath = getpref('Clintrajan','PythonPath');
        if ~isfile(pypath)
            fprintf('\nStored Python path no longer exists: %s\n', pypath);
            force_flag = true;
        end
    end

    % User input to select python installation path
    if isempty(pypath) || force_flag

        newpath = true;

        % Find Python installations (local function)
        python_paths = findPythonInstallationsInternal();

        if isempty(python_paths)
            fprintf('\nNo Python installations found automatically.\n');
            custom_path = input('Please enter your Python executable path manually: ', 's');
            if ~isempty(custom_path) && exist(custom_path, 'file')
                pypath = custom_path;
            else
                error('Valid Python path is required to continue.');
            end
        else
            % Display found Python installations
            n = numel(python_paths);
            fprintf('\nFound %d Python installation(s):\n',n);
            for i = 1:n
                fprintf('[%d] %s\n', i, python_paths{i});
            end
            
            % Let user choose which Python to use
            choice = input(sprintf('\nSelect Python installation [1-%d]: ', length(python_paths)));
            if isempty(choice) || choice < 1 || choice > n
                error('Invalid selection.');
            end
            pypath = python_paths{choice};
        end

    end
    
    if newpath
        % Check required packages
        flag = checkRequiredPackages(pypath);
        
        % SAVE THE PYTHON PATH FOR FUTURE USE
        setpref('Clintrajan', 'PythonPath', pypath);
        fprintf('\nPython path saved: %s\n', pypath);
    
        if flag
            fprintf('')
        else
            fprintf('You can now run the Clintrajan MATLAB functions with your selected Python.\n');
        end
    end
end

function python_paths = findPythonInstallationsInternal()
    % Internal function to find Python installations
    python_paths = {};
    
    % Check by OS type
    if ispc  % Windows
        % Try system Python first
        [status, ~] = system('where python');
        if status == 0
            [~, result] = system('where python');
            paths = strsplit(strtrim(result), '\n');
            for i = 1:length(paths)
                if ~isempty(paths{i}) && exist(paths{i}, 'file')
                    python_paths{end+1} = strtrim(paths{i});
                end
            end
        end
        
        % Check common locations
        [~,userpath] = system('echo %USERPROFILE%'); userpath = strtrim(userpath);
        potential_paths = {
            'C:\Python*\python.exe';
            'C:\Program Files\Python*\python.exe';
            'C:\ProgramData\Anaconda*\python.exe';
            fullfile(userpath,'AppData\Local\Programs\Python\Python*\python.exe');
            fullfile(userpath,'AppData\Local\anaconda*\python.exe');
            fullfile(userpath,'AppData\Local\miniconda*\python.exe');
            fullfile(userpath,'Anaconda*\python.exe');
            fullfile(userpath,'.conda\envs\*\python.exe')
        };
        
        for i = 1:length(potential_paths)
            files = dir(potential_paths{i});
            for j = 1:length(files)
                if ~files(j).isdir
                    full_path = fullfile(files(j).folder, files(j).name);
                    python_paths{end+1} = full_path;
                end
            end
        end
    else  % macOS and Linux
        % Try system Python first using 'which'
        [status, ~] = system('which python');
        if status == 0
            [~, result] = system('which python');
            python_paths{end+1} = strtrim(result);
        end
        
        [status, ~] = system('which python3');
        if status == 0
            [~, result] = system('which python3');
            python_paths{end+1} = strtrim(result);
        end
        
        % Check Homebrew and Anaconda common locations on macOS
        if ismac
            potential_paths = {
                '/usr/bin/python*';
                '/usr/local/bin/python*';
                '/opt/homebrew/bin/python*';
                '/opt/anaconda*/bin/python';
                '/opt/miniconda*/bin/python';
                fullfile(getenv('HOME'), 'anaconda*/bin/python');
                fullfile(getenv('HOME'), 'miniconda*/bin/python')
            };
            
            for i = 1:length(potential_paths)
                files = dir(potential_paths{i});
                for j = 1:length(files)
                    if ~files(j).isdir
                        full_path = fullfile(files(j).folder, files(j).name);
                        python_paths{end+1} = full_path;
                    end
                end
            end
        end
    end
    
    % Remove duplicates
    python_paths = unique(python_paths);
    
    % Verify paths actually work
    valid_paths = {};
    for i = 1:length(python_paths)
        if exist(python_paths{i}, 'file')
            cmd = sprintf('"%s" --version', python_paths{i});
            [status, ~] = system(cmd);
            if status == 0
                valid_paths{end+1} = python_paths{i};
            end
        end
    end
    
    python_paths = valid_paths;
end

function missing_flag = checkRequiredPackages(python_path)

    % Check for required Python packages
    required_packages = {'numpy', 'scipy', 'igraph', 'networkx', 'elpigraph',        'sklearn'};
    install_name =      {'numpy', 'scipy', 'igraph', 'networkx', 'elpigraph-python', 'scikit-learn'};
    np = numel(required_packages);
    missing_flag = false(1,np);
    
    fprintf('\nChecking for required Python packages...\n');
    
    for i = 1:np
        fprintf('  - %s : ',required_packages{i});

        cmd = sprintf('"%s" -c "import %s; print(''OK'')"', python_path, required_packages{i});
        [status, result] = system(cmd);
        
        if status ~= 0 || ~contains(result, 'OK')
            fprintf('missing\n');
            missing_flag(i) = true;
        else
            fprintf('found\n');
        end
    end
    
    for i = 1:np
        if missing_flag(i)
            if ispc
                missing_flag(i) = system(sprintf('"%s" -m pip install %s',python_path,install_name{i}));
            end
        end
    end

    missing_flag = any(missing_flag);

end