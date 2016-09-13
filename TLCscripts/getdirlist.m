% --------------------------
% Local function getdirlist
% --------------------------
function d=getdirlist(filefilter,fileExt)
 dirall = dir;
 s=size(dirall); s=s(1);
 clear d;count=0;
 for k=1:s
   if dirall(k).isdir
     % let us notice all directories
  %   count=count+1;
  %   d(count)=dirall(k);
  %   d(count).name=['<',d(count).name,'>'];
   elseif strcmp('.* ','.* ')
    a =dirall(k).name;
    if isempty(filefilter)
      if  ~isempty(strfind(lower(dirall(k).name), ['.',fileExt]))
        count=count+1;
        d(count)=dirall(k);
      end 
    elseif  ~isempty(strfind(lower(dirall(k).name), filefilter)) &&  ~isempty(strfind(lower(dirall(k).name), ['.',fileExt]))
     count=count+1;
     d(count)=dirall(k);
    end ;
   elseif  findstr([lower(dirall(k).name),' '],filefilter)>0
 %    count=count+1;
 %    d(count)=dirall(k);
   end;
 end;
% d(1).name='Q U I T';
return;
