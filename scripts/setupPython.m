function setupPython

% Check for version of Python
pe = pyenv;
if pe==0
    fprintf('Must install compatible version of Python!\n');
    return;
end

% Install required modules
mods = {'numpy','scipy','ftools','keras','tensorflow','scikit-image',...
    'SimpleITK','pynrrd'};
system(sprintf('pip install %s; ',mods{:}))

% Check for pycode path in Python search path
pypath = which('cmi');
pypath = fullfile(fileparts(pypath),'pycode');
if count(py.sys.path,pypath) == 0
    insert(py.sys.path,int32(0),pypath);
end