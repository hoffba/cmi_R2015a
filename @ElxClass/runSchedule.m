% ElxClass function
% Run Elastix Using Current Schedule
function stat = runSchedule(self,varargin)
stat = false;
if ~isempty(self.Schedule) && ...
        all(cellfun(@(x)exist(x,'file'),self.ImgInfo.FixedImage)) && ...
        all(cellfun(@(x)exist(x,'file'),self.ImgInfo.MovingImage))
    ns = length(self.Schedule);
    
    % First make sure all parameter files are saved
    for i = 1:ns
        if ~exist(self.Schedule(i).Name,'file')
            
        end
    end
    
    % Generate call to Elastix
    nf = length(self.ImgInfo.FixedImage);
    nfstr = {[]};
    if nf>1
        nfstr = num2cell(1:nf);
    end
    nfstr = [nfstr;self.ImgInfo.FixedImage];
    nm = length(self.ImgInfo.MovingImage);
    nmstr = {[]};
    if nm>1
        nmstr = num2cell(1:nm);
    end
    nmstr = [nmstr;self.ImgInfo.MovingImage];
    estr = ['/opt/elastix/bin/elastix',...
            sprintf(' -f%u %s',nfstr{:}),...
            sprintf(' -m%u %s',nmstr{:}),...
            sprintf(' -p %s',self.Schedule(:).Name),...
            ' -out ',self.wdir];
    if ~isempty(self.ImgInfo.FixedMask)
        estr = [estr,' -fMask ',self.ImgInfo.FixedMask];
    end
    if ~isempty(self.ImgInfo.MovingMask)
        estr = [estr,' -mMask ',self.ImgInfo.MovingMask];
    end
    [stat,str] = system(estr);
end
