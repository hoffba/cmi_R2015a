function [h,fv] = surfmap(varargin)
% h = surfmap(BW,voxsz)
% h = surfmap(BW,voxsz,V,clim)
% h = surfmap(fv,V,clim)
% h = surfmap(fv,clim)

V = [];
switch nargin
    case 2
        if isstruct(varargin{1})
            fv = varargin{1};
            clim = varargin{2};
            if isfield(fv,'facevertexcdata')
                V = fv.facevertexcdata;
            end
        elseif islogical(varargin{1})
            fv = mask2surf(varargin{:});
        else
            error('Invalid input.')
        end
    case 3
        fv = varargin{1};
        V = varargin{2};
        clim = varargin{3};
    case 4
        fv = mask2surf(varargin{1:3});
        clim = varargin{4};
        if isfield(fv,'facevertexcdata')
            V = fv.facevertexcdata;
        end
    otherwise
        error('Invalid inputs.');
end
opts = {'EdgeColor','none',...
        'FaceLighting','gouraud',...
        'AmbientStrength',0.5,...
        'DiffuseStrength',0.5,...
        'SpecularStrength',0.3,...
        'SpecularExponent',50,...
        'BackFaceLighting','reverselit'};
nv = length(V);
if nv == 0
    opts = [{'Facecolor','flat'},opts];
elseif nv == length(fv.faces)
    opts = [{'FaceVertexCData',V,'Facecolor','flat'},opts];
elseif nv == length(fv.vertices)
    opts = [{'FaceVertexCData',V,'Facecolor','interp'},opts];
else
    error('Invalid CData input length.');
end

hf = figure('Colormap',jet(128),'Position',[500 500 1200 800],'Units','normalized');
ha1 = axes(hf,'Position',[0 0 0.5 1],'CLim',clim);
la1 = lightangle(60,-30);
ha2 = axes(hf,'Position',[0.5 0 0.5 1],'CLim',clim);
la2 = lightangle(240,-30);

p1 = patch(ha1,'Vertices',fv.vertices,'Faces',fv.faces,opts{:});
axis(ha1,'equal','off','tight');
view(ha1,90,0);

p2 = patch(ha2,'Vertices',fv.vertices,'Faces',fv.faces,opts{:});
axis(ha2,'equal','off','tight');
view(ha2,270,0);

h = struct('hfig',hf,'axes1',ha1,'axes2',ha2,...
    'lightAngle1',la1,'lightAngle2',la2,'patch1',p1,'patch2',p2);


