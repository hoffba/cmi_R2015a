function bvalue = getbvalue(info)
% this function is used to retrieve bvalue (Philips)
% TLChenevert, UMICH    Feb 2006
% TLChenevert, UMICh Feb 2009:  Update to read bvalue from GE and SIEMENS dicom.
% TLC 7/22/2011 GE bvalues on new platform may have 1e9 DC added for some reason.
% TLC 8/15/2011.  Temporarily hardcode manufacturer = 'GE' to read De-ID HFH GE data Search by "8/15/2011".
% TLC 5/30/2012.  Changes to handle uint8 GE bvalues.
% TLC 7/07/2014.    Add Agilent as a mfg
% MD 20141006: uint8 may also need -1e9 (moved out of "esle" )
% TLC 20150827.  Switch to tryGetField(info, 'Manufacturer','UNK')
% TLC 20150904.  Add manufacturer = Imaging Biometrics LLC.
% MD 20151009: Edited to use public tag for Philips "DiffusionBValue" 
% MD 20151012: get Siemens b-val from "SequenceName" (0018;0024) per D.Newitt suggestion
% TLC 20151226: Siemens ADC maps have info.SequenceName = '*ep_b0_1000';
% which will be returned as '[]'. This causes problems in sortedre.
% Therefore do final check against this (should apply to any vendor).
% TLC 20160629:  Resort to private Philips bvalue tag if public tags not found.

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
% % manufacturer
% % Temp fight thru old GE data ...
% % manufacturer = 'GE'; % 8/15/2011.
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
        % bvalue may not exist  
        bvalexist = isfield(info,'DiffusionBValue'); % 1st try public tag
        bvalexistalt = isfield(info,'Private_2001_1003'); % 2nd try private tag
       % bvalexist = isfield(info,'Private_2001_1003');
       % 20160629 start.  First, blot out old ...
%         if (bvalexist == 0)
%             
%             bvalue = 0;
%         else
%             bvalue = info.DiffusionBValue; %MD 20151009
%             %bvalue = info.Private_2001_1003;    
%             if (isa(bvalue,'uint8') == 1)
%                     bvalue = hex2float(arr2str(bvalue));
%             end % if isa.  Must already be a 'double'.           
%         end % if bvalexist
        if (bvalexist == 1)
            bvalue = info.DiffusionBValue;
            if (isa(bvalue,'uint8') == 1)
                    bvalue = hex2float(arr2str(bvalue));
            end % if isa.  Must already be a 'double'.           
        elseif (bvalexistalt == 1)
            bvalue = info.Private_2001_1003;
            if (isa(bvalue,'uint8') == 1)
                bvalue = hex2float(arr2str(bvalue));
            end % if isa.  Must already be a 'double'.
        else
            bvalue = 0;
        end % if bvalexist
            
    
    case 2 % Siemens
        % bvalue may not exist    
        bvalexist = isfield(info,'Private_0019_100c');
        if (bvalexist == 0)
            seqnm = tryGetField(info, 'SequenceName','UNK');
            kb = strfind(seqnm, '_b'); % assuming standrad name: "*epi_b###t"
            if ( ~isempty(kb) && ((kb+2) < (length(seqnm)-1)) )
               bvalue = str2num(seqnm((kb+2):(end-1))); % getting b-val from "*epi_b###t"
            else
                bvalue = 0;
            end
        else
            bvalue = info.Private_0019_100c; % Assume double already, check later if not.
            if (isa(bvalue,'uint8') == 1)
                bvalue = str2num(char(bvalue'));
            end % if isa          
        end % if bvalexist
        
    case 3 % GE     
%         bvalue = 0; % Find out later where GE bvalues are
        bvalexist = isfield(info,'Private_0043_1039');
        if (bvalexist == 0)
            bvalue = 0;
        else
            zz = info.Private_0043_1039;
            if (isa(zz,'uint8') == 1)
                % disp('bvalues coded as uint8 characters')
                str = char(zz'); % UAB DWI via TRIAD/ACRIN has bvalues in strings, eg. "800\8\0\0" for bvalue = 800.
                istr = strfind(str,'\');
                bvalue = str2num(str(1:(istr(1) - 1)));      
            else
                bvalue = zz(1,1);
            end % if isa uint8
            %%% MD20141006: uint8 may alo need -1e9
                if (bvalue > 1e6)
                    %disp('***');
                    bvalue = bvalue - 1e9;
                end % if > 1e6
        end % if bvalexist    
    otherwise
        bvalue = 0;

end % end switch mfg

% 20151226 start ... Last check against returning empty
if( isempty(bvalue) )
    bvalue = 0;
end % if isempty
% ... end 20151226

return;

function tt2=arr2str(arr)
tt=dec2hex(arr,2);
tt2=[tt(4,1) tt(4,2) tt(3,1) tt(3,2) tt(2,1) tt(2,2) tt(1,1) tt(1,2)];
return;

function x=hex2float(s)
% Converts a hexadecimal 8-byte character string into the equivalent 32-bit floating
% point number. The input can be a cell array of strings, the strings are
% padded to the right with '0' if need be.
% The rules obey the IEEE-754 standard, although +/- 0 is not supported. 
% The first 9 lines of code are taken from 'hex2num', with slight modification
% B.Babic 29 Nov 1999
if iscellstr(s), s = char(s); end;
if ~isstr(s)
error('Input to hex2float must be a string.');
end;
if isempty(s), x = []; return, end;

[row,col] = size(s);
blanks = find(s==' '); 
% zero pad the blanks
if blanks, s(blanks) = '0'; end ;

d=48*ones(row,8);
i=find(s);
d(i)=abs(lower(s));
if ~(all((d<58&d>47)|(d<103&d>96)))
error('Input should only contain 0-9, a-f, A-F');
end;

d=(d-48).*(d<58)+(d-87).*(d>58);
s=(-1).^(d(:,1)>7);
e=32*d(:,1)-256*(d(:,1)>7)+2*d(:,2)+(d(:,3)>7);
f=[2*d(:,3:7)-16*(d(:,3:7)>7)+(d(:,4:8)>7), 2*d(:,8)-16*(d(:,8)>7)];
% f=f./16.^[1:6]
f=f*[6.250000000000000e-002;3.906250000000000e-003;2.441406250000000e-004;...
1.525878906250000e-005;9.536743164062500e-007;5.960464477539063e-008];

i1=find(e<255&e>0); % valid number
i2=find(e==255&f); % NaN
i3=find(e==255&~f); % Inf's
i4=find(e==0); % denormalized

x=zeros(1,row);
if i1,x(i1)=pow2(1+f(i1),e(i1)-127).*s(i1);end;
if i2,x(i2)=nan;end;
if i3,x(i3)=s(i3).*inf;end;
if i4,x(i4)=pow2(f(i4),e(i4)-126).*s(i4);end;
return;
