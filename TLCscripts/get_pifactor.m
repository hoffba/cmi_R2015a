function [pi_used pi_fctr]  = get_pifactor(ii, vndr)
% 
% Auxiliary function to get parallel imaging info, depending on "vendor-code" input string.
% 
% MD 20160523: extracted from "build_acringroup" and modified with "tryGetField (for safety)"
%
%vndr  %debug
% Parallel Imaging Section:
if (strcmp(vndr,'PMS'))

    if( isfield(ii,'ParallelAcquisition') ) % Good, parallel imaging info is at top level - use it.
        pi_used = ii.ParallelAcquisition; % "YES" or "NO"
%        if( isfield(ii,'ParallelReductionFactorInPlane') )
        pi_fctr = tryGetField(ii, 'ParallelReductionFactorInPlane', 1);
%        else
%           pi_fctr = 1;
%        end % isfield 'ParallelReductionFactorInPlane'
    elseif( isfield(ii,'Private_2005_140f')) % 20140312
        if (isstruct(ii.Private_2005_140f)) % MD 20160523 : check for proper format
            pi_used = tryGetField(ii.Private_2005_140f.Item_1, 'ParallelAcquisition', 'UNK'); % New for ACRIN6698
            if (pi_used(1:2) == 'YE')
                pi_fctr = tryGetField(ii.Private_2005_140f.Item_1, 'ParallelReductionFactorInPlane', 1); % New for ACRIN6698
            else % Then is "NO"
                pi_fctr = 1;
            end % if Group.parallel 
        else % not a structure
            pi_used = 'UNK';
            pi_fctr = 0; 
        end
    else
        pi_used = 'UNK';
        pi_fctr = 0;
    end % isfield 'ParallelAcquisition'
    % ... End 20140822
    
elseif (strcmp(vndr,'GMS'))
    %Private_0043_1083 is GE ASSET scan reductionfactor
    if (isfield(ii,'Private_0043_1083'))
        zz = ii.Private_0043_1083;
        if (isa(zz,'uint8') == 1) % Start MAY31_2012
            % disp('ASSET reduction factor coded as uint8 characters')
            % UAB DWI via TRIAD/ACRIN has reduction in strings, eg. "0.5\1" for reductionfactor = 0.5 or ASSET = 2.
            str = char(zz');
            istr = strfind(str,'\');
            reductionfactor = str2num(str(1:(istr(1) - 1)));
        else
            reductionfactor = zz;
        end % if isa uint8   % End MAY31_2012     
        if ( (reductionfactor(1) < 1) && (reductionfactor(1) ~= 0) ) 
            pi_used = 'YES';
            pi_fctr = 1/reductionfactor(1);
        else
            pi_used = 'NO';
            pi_fctr = 1;
        end % if reductionfactor
    else
        pi_used = 'UNK';
        pi_fctr = 0;
    end % if isfield

% Start 20121210 ...>>>>>
elseif (strcmp(vndr,'SMS')) % Added Siemens Parallel Img info on APR29_2012
    %isfield(ii,'Private_0051_1011') % debug
    if (isfield(ii,'Private_0051_1011')) % I think this is only here in Siemens when parallel img applied
        junk1 = ii.Private_0051_1011; % This is a char array, like 'p2'
        if (isnumeric(junk1)) junk1 = char(junk1); end;
        szj1 = size(junk1);
        % junk1
        if (szj1(2) > szj1(1))
            junk2 = junk1(1,2); % 2nd character looks like the GRAPPA factor
        elseif (szj1(1) > szj1(2))
            junk2 = junk1(2,1); % 2nd character looks like the GRAPPA factor
        else
            disp('WARNING: no GRAPPA provided');
            junk2 = junk1(1,1);
        end
        junk3 = str2num(junk2);

        if ( (junk3 >= 1) ) 
            pi_used = 'YES';
            pi_fctr = junk3;
        else
            pi_used = 'NO';
            pi_fctr = 1;
        end % if reductionfactor
    else
        pi_used = 'UNK';
        pi_fctr = 0;
    end % if isfield End APR29_2012 insert 
else % Need to look into Siemens later
    pi_used = 'UNK';
    pi_fctr = 0;
end % if vendor-code

return