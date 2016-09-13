function mfg = getmfg(info)
% This function is used to retrieve manufacturer that subsequenlty drives
% dicom tag selection for other image properties.  This script serves
% getbvalue, getdiffgrad, and getimagetype.
% 20150909  TLChenevert
% 20160520  TLC Deal with condition that field exists, but is empty
% 20160523  TLC NRG case 204 Varian mixed in with Philips data

% Initialize all manufacturer flags to 0:
isPhilips = 0;
isSIEMENS = 0;
isGE = 0;
isAgilent = 0;
isVarian  = 0; % 20160523
isBiometrics = 0; % 20150904
isUnknown = 1; % 20150827

% Read manufacturer:
% manufacturer = getdicom_field_str(info, 'Manufacturer'); %8/15/2011.
manufacturer = tryGetField(info, 'Manufacturer','UNK'); % 20150827
%manufacturer = 'GE'; % 8/15/2011.

if (isempty(findstr(manufacturer ,'Philips')) == 0)
    % Then is Philips data
    isPhilips = 1;
    mfg = 1;
end
if (isempty(findstr(manufacturer ,'SIEMENS')) == 0)
    % Then is SEIMENS data
    isSIEMENS = 1;
    mfg = 2;
end
if (isempty(findstr(manufacturer ,'GE')) == 0)
    % Then is GE data
    isGE = 1;
    mfg = 3;
end
if (isempty(findstr(manufacturer ,'Agilent')) == 0)
    % Then is Agilent/Varian data
    isAgilent = 1;
    mfg = 10; % Skip over other future human scanners
end
% 201500904 start ...
if (isempty(findstr(manufacturer ,'Biometrics')) == 0)
    % Then is unknown data
    isBiometrics = 1;
    mfg = 11; % Skip over other future human scanners
end
% 20160523 start ...
if (isempty(findstr(manufacturer ,'Varian')) == 0)
    % Then is Agilent/Varian data
    isVarian = 1;
    mfg = 12; % Rad-Onc system
end
% ... 20160523 end.

% ... 20150904 end.
% % 20150827 start ...
% if (isempty(findstr(manufacturer ,'UNK')) == 0)
%     % Then is unknown data
%     isUnknown = 1;
%     mfg = 99; % Skip over other future human scanners
% end
% % ... 20150827 end.

% 20160520 start ...
if ( (isempty(findstr(manufacturer ,'UNK')) == 0) || (isempty(manufacturer)) )
    % Then is unknown data
    isUnknown = 1;
    mfg = 99; % Skip over other future human scanners
end
% ... 20160520 end.

return;
