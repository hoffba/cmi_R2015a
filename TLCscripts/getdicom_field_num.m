function value = getdicom_field_num(info, field)
% this function is used to retrieve dicom field
% output augument is returned as a number
% HG    11/23/05 
   
    % check if the field exists or not    
    fieldexist = isfield(info, field);
    if (fieldexist == 0)
        value = 0;  %default
    else
        value = info.(field);
    end

return;