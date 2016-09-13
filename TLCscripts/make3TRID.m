function [def3TRID] = make3TRID(dicinfo,diid)
% Create 3T Research ID code for given series based on dicom info.
% diid = 'di' (default) Use de-id transformation and have "di" as 1st characters.
% diid = 'id' Do not transform and have "id" as 1st characters.
% 20150627 TLChenevert

if (nargin < 2)
    diid = 'di';
end

if (diid == 'di')
    xform = 1;
elseif (diid == 'id')
    xform = 0;
else
    disp('Aborting, Second argument should be either "di" or "id". ');
    return;
end

sdate = tryGetField(dicinfo,'SeriesDate','00000000'); % 20150204 returns "00000000" if no field present
fndum1 = str2num(sdate(3:4)) + 3*xform;
if (fndum1 < 10)
   fndum11 = ['0' num2str(fndum1)];
else
   fndum11 = num2str(fndum1);
end
fndum2 = str2num(sdate(5:6)) + 2*xform;
if (fndum2 < 10)
   fndum22 = ['0' num2str(fndum2)];
else
   fndum22 = num2str(fndum2);
end
fndum3 = str2num(sdate(7:8)) + 1*xform;
if (fndum3 < 10)
   fndum33 = ['0' num2str(fndum3)];
else
   fndum33 = num2str(fndum3);
end
% fndum4 = dicinfo.SeriesTime(1:4);
fndum4 = tryGetField(dicinfo,'SeriesTime','0000'); % 20150204 returns "0000" if no field present
fndum4 = fndum4(1:4); % 20150204
fndum5 = dicinfo.SeriesNumber;
%fndum6 = [fndum11 fndum22 fndum33 num2str(fndum4) 's' num2str(fndum5) '_0001'];
fndum6 = [diid fndum11 fndum22 fndum33 num2str(fndum4) 's' num2str(fndum5)];
def3TRID = fndum6;

end

