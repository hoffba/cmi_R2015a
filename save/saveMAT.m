function status = saveMAT(fname,img,label,fov)
%% Save .mat file with fields of: "img", "label", and "fov"
status = 0;
if (nargin==4) && ~isempty(img) && ~isempty(label) && ~isempty(fov)
    label = {strtok(label,'label= ""')};
    save(fname,'img','label','fov');
    status = 1;
end