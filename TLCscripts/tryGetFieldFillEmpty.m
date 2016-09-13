function val = tryGetFieldFillEmpty(s, field, dftVal)
% TLC 20150204.  Externalized from "build_nii.m"  Originally from Xiangrui Li available
% on MatLab User Community website
% 20150716 DM TC.  Make sure returned "val" is not a structure, but is a field.
% 20151106 TLC.  Variation of "tryGetField" that fills empty fields with
% dftVal.  Unlike tryGetField, the 3rd argument, dftVal, must be provided.

if isfield(s, field)
    val = s.(field);
    
    % 20150716 Start ...
    if (isstruct(val))
            fn = fieldnames(val);
            vv = getfield(val,char(fn(1)));
            val = vv;
    end
    % ... End 20150716
    if ( isempty(val) )
        val = dftVal;
    end
else
    val = dftVal;
end


