function ppd = getppd(info)
% this function is used to retrieve prepause delay (Philips)
% TLChenevert, UMICH    Feb 2006

    % prepause may not exist    
    ppdexist = isfield(info,'Private_2001_101b');
    if (ppdexist == 0)
        ppd = 0;
    else
        ppd = info.Private_2001_101b;
        if (isa(ppd,'uint8') == 1)
            ppd = hex2float(arr2str(ppd));
        end % if isa.  Must already be a 'double'.
    end

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
