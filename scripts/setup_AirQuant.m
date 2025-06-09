function setup_AirQuant
% Installs PTK and AirQuant in directories inside the CMI program

install_path = fileparts(which('cmi'));
cmdstr = ['cd "',install_path,'"',...
    ' & git clone https://github.com/ashkanpakzad/AirQuant.git',...
    ' & git clone https://github.com/ashkanpakzad/pulmonarytoolkit.git',...
    ' & cd pulmonarytoolkit',...
    ' & git checkout nifti-gz-load-import'];
[stat,outstr] = system(cmdstr);
if stat

else
    error(['Failed to set upt AirQuant : ',outstr]);
end