function imty = getimagetype(info)
% this function is used to retrieve bvalue (Philips)
% TLChenevert, UMICH    Feb 2006
% TLChenevert, UMICh Feb 2009:  Update to read bvalue from GE and SIEMENS dicom.
% TLC 7/22/2011 GE bvalues on new platform may have 1e9 DC added for some reason.
% TLC 8/15/2011.  Temporarily hardcode manufacturer = 'GE' to read De-ID HFH GE data Search by "8/15/2011".
% TLC 5/30/2012.  Changes to handle uint8 GE bvalues.
% TLC 20140221: Need image type (egs, Mag, Real, Phase, ADC, FA, Velocity, ... for proper 4th dimension sorting.
%               Returns a numerical value associated with image type:
%                   0  = unknown
%                   1  = magnitude
%                   2  = real
%                   3  = imaginary
%                   4  = phase
%                   5  = adc
%                   6  = fa
%                   7  = velocity
%                   8  = water
%                   9  = in-phase
%                   10 = out-phase
%                   11 = fat
%                   12 = T1
%                   13 = T2
%                   14 = OTHER
%                   15 = SCREEN SHOT
%                   16 = SECONDARY
%                   17 = FF         Philips Fat Fraction
%                   18 = T2_STAR    Philips T2*
%                   ... add more as needed
%                   99 = currently anthing else
% TLC 20140707: Add Agilent
% TLC 20141223: Add logic for "1" backslash and type = SECONDARY
% TLC 20150303: Add cases for Philips QuantDixon (FatFraction and T2_Star)
% TLC 20150827.  Switch to tryGetField(info, 'Manufacturer','UNK')
% TLC 20150904.  Add manufacturer = Imaging Biometrics LLC.
% TLC 20151218.  Needed work-around for Siemens "ImageType", when they add 4th backslash.

% % Initialize all manufacturer flags to 0:
% isPhilips = 0;
% isSIEMENS = 0;
% isGE = 0;
% isAgilent = 0;
% isBiometrics = 0; % 20150904
% isUnknown = 1; % 20150827
% 
% % Read manufacturer:
% % manufacturer = getdicom_field_str(info, 'Manufacturer'); %8/15/2011.
% manufacturer = tryGetField(info, 'Manufacturer','UNK'); % 20150827
% %manufacturer = 'GE'; % 8/15/2011.
% 
% if (isempty(findstr(manufacturer ,'Philips')) == 0)
%     % Then is Philips data
%     isPhilips = 1;
%     mfg = 1;
% end
% if (isempty(findstr(manufacturer ,'SIEMENS')) == 0)
%     % Then is SEIMENS data
%     isSIEMENS = 1;
%     mfg = 2;
% end
% if (isempty(findstr(manufacturer ,'GE')) == 0)
%     % Then is GE data
%     isGE = 1;
%     mfg = 3;
% end
% if (isempty(findstr(manufacturer ,'Agilent')) == 0)
%     % Then is Agilent/Varian data
%     isAgilent = 1;
%     mfg = 10; % Skip over other future human scanners
% end
% % 201500904 start ...
% if (isempty(findstr(manufacturer ,'Biometrics')) == 0)
%     % Then is unknown data
%     isBiometrics = 1;
%     mfg = 11; % Skip over other future human scanners
% end
% % ... 20150904 end.
% % 20150827 start ...
% if (isempty(findstr(manufacturer ,'UNK')) == 0)
%     % Then is unknown data
%     isUnknown = 1;
%     mfg = 99; % Skip over other future human scanners
% end
% % ... 20150827 end.

mfg = getmfg(info);

switch mfg
    case 1 % Philips
        % ImageType
        if (isfield(info,'ImageType'))
            imgtyp = info.ImageType;
        else
            imgtyp = 'UNK\UNK\UNK\UNK';
        end
         % or could use: imgtyp = tryGetField(info, 'ImageType', 'UNK\UNK\UNK\UNK')
        
    
    case 2 % Siemens
        % ImageType
        if (isfield(info,'ImageType'))
            imgtyp = info.ImageType;
            % 20151218 start ...
            bs = findstr(imgtyp,'\'); % count backslashes.
            [junk, ln2] = size(bs);
            [junk, ln3] = size(imgtyp);
            % Need to limit to 3 backslashes
            if (numel(bs) > 3) % Check that 3rd and 4th backslash exists.
                junk = imgtyp(1:(bs(4)-1)); % Clip off at last bs
                imgtyp = junk;
            end
            % ... 20151218 end

        else
            imgtyp = 'UNK\UNK\UNK\UNK';
        end
        
    case 3 % GE     
        % ImageType
        if (isfield(info,'ImageType'))
            imgtyp = info.ImageType;
        else
            imgtyp = 'UNK\UNK\UNK\UNK';
        end
        
    otherwise
        if (isfield(info,'ImageType'))
            imgtyp = info.ImageType;
        else
            imgtyp = 'UNK\UNK\UNK\UNK';
        end

end % end switch mfg


bs = findstr(imgtyp,'\'); % count backslashes.

[junk, ln2] = size(bs);
[junk, ln3] = size(imgtyp);


if (numel(bs) > 3) % Check that 3rd and 4th backslash exists.
    imagetype = imgtyp((bs(3)+1):(bs(4)-1)); % expect this is character strings.
elseif (numel(bs) > 2) % Have GE cases with only 3
    imagetype = imgtyp((bs(2)+1):(bs(3)-1)); % 'OTHER' and 'SCREEN SAVE'
elseif (numel(bs) == 2) % Have GE cases with only 2
    imagetype = imgtyp((bs(2)+1):end); % 'OTHER' and 'SCREEN SAVE'
elseif ((numel(bs) == 1)) % Have Philips cases with only 1
    imagetype = imgtyp((bs(1)+1):end); % Eg Philips 'ORIGINAL\SECONDARY' becomes 'SECONDARY'
else
    imagetype = 'UNK';
end % if numel

switch imagetype % Most of these are based on Philips, add more as discovered.  Some strings are assumed - fix later when correct string found.
    case 'UNK'          % Default in event "ImageType" tag is not present
        imty = 0; 
    case 'M'            % Magnitude image
        imty = 1; 
    case 'R'            % Real image
        imty = 2;
    case 'I'            % Imaginary image
        imty = 3;
    case 'P'            % Phase image
        imty = 4;
    case 'ADC'          % ADC map - caution, do not yet know units, but likely "1000" = 1x10-3mm2/s
        imty = 5;
    case 'FA'           % Fractional Anisotropy, do not know units but likely FA=1 displayed as 1000.
        imty = 6;
    case 'V'            % Velocity image - units likely mm/s or cm/sec x 1000
        imty = 7;
    case 'W'            % Water-only image, from Dixon recon
        imty = 8;
    case 'IP'           % In-Phase image, from Dixon recon
        imty = 9;
    case 'OP'           % Out-of-Phase image, from Dixon recon
        imty = 10;
    case 'F'            % Fat-only image, from Dixon recon
        imty = 11;
    case 'T1'           % T1 map, units unknown but likely 1000 = 1000ms
        imty = 12;
    case 'T2'           % T2 map, units unknown but likely 100 = 100ms
        imty = 13;
    case 'OTHER'        % Have GE case where "OTHER" used - no meaning here.
        imty = 14;
    case 'SCREEN SHOT'  % Have GE case where secondary capture screen shot used.
        imty = 15;
    case 'SECONDARY'    % Have Philips case where secondary capture used for scanned docs.
        imty = 16;
    case 'FF'           % Philips Fat Fraction from Quantitative mDixon. *** Read as 'dv' in percent.  Pixels beyond FOV = -100 and artifact=100-200***.  
        imty = 17;
    case 'T2_STAR'      % Philips T2* from Quantitative mDixon.  *** Read as 'dv' in ms ***.  
        imty = 18;
    % Add more cases as discovered....    
    otherwise
        imty = 99;
end

return;
