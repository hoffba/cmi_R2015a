% CMI script
function AddJobMouseLung(cmiObj)
% Starts batch processes for performing NOVA filter and saving as MHD
%   Input: cmiObj = CMIclass object containing current settings

% Check to see if Parallel Processing is available
if isempty(license('inuse','Distrib_Computing_Toolbox'))
    TF = license('checkout','Distrib_Computing_Toolbox');
else
    TF = true;
end
if TF
    
    % Set image labels:
    answer = inputdlg('Which image is the end-expiratory?','Exp Image',1,{'1'});
    if isempty(answer)
        return;
    else
        nv = cmiObj.img.dims(4);
        t0 = round(str2double(answer));
        if isnan(t0) || (t0<1) || (t0>cmiObj.img.dims(4))
            error(['Index out of available range: ',num2str(t0)]);
        else
            cmiObj.img.setLabel(1:nv,...
                circshift(cellfun(@(x)sprintf('%02u',x),num2cell(0:nv-1),...
                    'UniformOutput',false),t0-1,2));
        end
    end
    
    % Initialize batch inputs:
    fov = cmiObj.img.voxsz.*cmiObj.img.dims(1:3);
    bname = cmiObj.img.name(1:end-12);
    svdir = cmiObj.img.dir; % Save processed images to same directory
    hw = waitbar(0,'Starting batch processes ...');
    
    % Close pool before you can start a batch process:
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        % Can't run batch with open matlabpool
        delete(poolobj);
    end
    
    % Create batch process for each image
    % * Images too large to do all at once
    for i = 1:nv
        fname = [bname,cmiObj.img.labels{i},'.mhd'];
        j = batch(@jobProcMouseLung,0,{cmiObj.img.mat(:,:,:,i),...
                    fov,cmiObj.img.labels(i),...
                    fname,svdir});
        disp([fname,' : Job # ',num2str(j.ID)]);
        waitbar(i/nv,hw);
    end
    delete(hw);
else
    disp('You do not have a Distrib_Computing_Toolbox license!');
end