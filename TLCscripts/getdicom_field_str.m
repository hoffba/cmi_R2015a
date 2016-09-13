function value = getdicom_field_str(info, field)
% this function is used to retrieve dicom field
% output augument is returned as a number
% HG    11/23/05 
% DM & TLC May2013 in anticipation of Matlab12
   
    % check if the field exists or not    
    fieldexist = isfield(info, field);
    if (fieldexist == 0)
        value = 'Not Available';  %default
    else
        value = info.(field);
        if (isstruct(value))
            fn = fieldnames(value);
            vv = getfield(value,char(fn(1)));
            value = vv;
        end

    end

return;