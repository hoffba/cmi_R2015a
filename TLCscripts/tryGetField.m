function val = tryGetField(s, field, dftVal)
% TLC 20150204.  Externalized from "build_nii.m"  Originally from Xiangrui Li available
% on MatLab User Community website
% 20150716 DM TC.  Make sure returned "val" is not a structure, but is a field.

if isfield(s, field)
    val = s.(field);
    
    % 20150716 Start ...
    if (isstruct(val))
            fn = fieldnames(val);
            vv = getfield(val,char(fn(1)));
            val = vv;
    end
    % ... End 20150716  
    
elseif nargin>2
    val = dftVal;
else
    val = [];
end


