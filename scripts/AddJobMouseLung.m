% CMI script
function AddJobMouseLung(cmiObj)
%   Input: cmiObj = CMIclass object containing current settings

% Check to see if Parallel Processing is available
if isempty(license('inuse','Distrib_Computing_Toolbox'))
    TF = license('checkout','Distrib_Computing_Toolbox');
else
    TF = true;
end
if TF
    if strfind(cmiObj.img.name,'-phase-')
        % Assumes fname is in format: ...-phase-00
        fname = cmiObj.img.name(1:end-9);
    else
        fname = cmiObj.img.name;
    end
    fov = cmiObj.img.voxsz.*cmiObj.img.dims(1:3);
    svdir = uigetdir(pwd,'Select directory to save DICOMs');
    if svdir~=0
        hw = waitbar(0,'Starting batch process ...');
        poolobj = gcp('nocreate');
        if ~isempty(poolobj)
            % Can't run batch with open matlabpool
            delete(poolobj);
        end
        j = batch(@jobProcMouseLung,1,{cmiObj.img.mat,fov,cmiObj.img.labels,fname,svdir});
        delete(hw);
        disp([fname,' : Job # ',num2str(j.ID)]);
    end
else
    disp('You do not have a Distrib_Computing_Toolbox license!');
end