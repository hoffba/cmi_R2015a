function assignedfieldvalue = getexistfield(structinput,possiblefield)
% Extract "possiblefield" value from structure "structinput" if it exists.
% If it does not exist, assing 'NA' value to "assignedfieldvalue".  Main use of
% this script is to extract values from dicominfo-derived structures.
% Often the desired field is not assinged in the structure and rather than
% having the script crash, it is preferred to return an 'NA' value.
% TLChenevert 8/23/2011.


    
fieldexist = isfield(structinput, possiblefield);
if (fieldexist == 0)
    assignedfieldvalue = 'NA';
else
    assignedfieldvalue = getfield(structinput,possiblefield);               
end % if fieldexist