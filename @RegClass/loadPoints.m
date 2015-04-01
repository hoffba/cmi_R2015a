% RegClass function
function loadPoints(self,i,fname)
% Load points from Elastix .txt file with format:
%       <"index" or "point"> (no quotes)
%       < number of points >
%       point1_x point1_y [point1_z]
%       point2_x point2_y [point2_z]
%           ...
% Inputs:
%       i - (1) Reference image, (2) Homologous image
%           default = selected radiobutton in imgReg GUI
%       fname - full file name of .txt file to save
%           default = user selects name/location
% Syntax:
%       fname = RegObj.loadPoints;
%       fname = RegObj.loadPoints(i);
%       fname = RegObj.loadPoints(i,fname);

if ((nargin==3) && ishandle(i)) || (nargin==1)
    i = find(cell2mat(get([self.h.radRef,self.h.radHom],'Value')),1);
end
p = self.ctrpts(i);
if ~isempty(p)
    if (nargin<3) || ~ischar(fname)
        [fname,fpath] = uiputfile('*.txt','Save Points');
        if fname~=0
            fname = fullfile(fpath,fname);
        else
            fname = '';
        end
    end
    if ~isempty(fname)
        fid = fopen(fname,'w');
        if fid>2
            fprintf(fid,'points\n%u\n',size(p,1));
            fprintf(fid,'%.8f %.8f %.8f\n',p');
            fclose(fid);
        end
    end
end


