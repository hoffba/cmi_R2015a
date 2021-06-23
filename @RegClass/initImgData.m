% RegClass function
function initImgData(self,varargin)

tObj = varargin{1};
i = (tObj.h.mainFig==self.cmiObj(2).h.mainFig)+1;
self.points{i} = [];
% Delete any exist point plots associated with image:
if ~isnan(self.hpts(i)) && ishandle(self.hpts(i))
    delete(self.hpts(i))
end

fvoxsz = self.cmiObj(1).img.voxsz([2,1,3]);
fdims = self.cmiObj(1).img.dims([2,1,3]);
if (i==1) % Reference image loaded
    self.elxObj.setTx0par('Size',fdims,'Spacing',fvoxsz);
else % If Homologous image was loaded, reset output directory based on name
    if nargin>1
        tdir = varargin{2};
    else
        tdir = fileparts(self.odir);
    end
    nname = self.cmiObj(2).img.name;
    ndir = self.cmiObj(2).img.dir;
    if ~isempty(nname)
        tname = nname;
    else
        tname = 'Unknown';
    end
    if ~isempty(ndir)
        tdir = ndir;
    end
    self.setOdir(fullfile(tdir,['elxreg_',tname,filesep]));
    
    % Initializes the transform
    self.elxObj.setTx0([1 0 0 0 1 0 0 0 1 0 0 0],fvoxsz,fdims);
    
    % Update GUI objects:
    self.showTx0;
    
end