% ImageClass function
function ok = imgMath(self,str)
% Input: str -- string of equation for calculation using ImageClass 
%               structure "self" as if calculating from within the object

ok = true;
if (nargin==1) || isempty(str) || ~ischar(str)
    str = inputdlg({'Type equation:','New Label'},'Image Math',[10,60;1,60],{'','Math'},'on');
    if isempty(str)
        ok = false;
    end
end
if ok
    try
        timg = eval(str{1});
        [d(1),d(2),d(3),d(4)] = size(timg);
    catch err
        assignin('base','err',err)
        disp(err.message);
        ok = false;
    end
    if ok && all(d(1:3)==self.dims(1:3))
        timg(isnan(timg)|isinf(timg)) = 0;
        self.mat(:,:,:,end+1) = timg;
        self.valExt(end+1,:) = [squeeze(min(min(min(timg,[],3),[],2),[],1)),...
                                squeeze(max(max(max(timg,[],3),[],2),[],1))];
        self.labels(end+1) = str(2);
        self.scaleM(end+1) = 1;
        self.scaleB(end+1) = 0;
        self.thresh(end+1,:) = [-100000,100000];
        self.dims(4) = size(self.mat,4);
    else
        ok = false;
    end
end