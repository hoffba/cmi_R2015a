function h = savesurf(vdm,flag,fname)
%,f,v,V,clim,label,logdisp,fname)

h = [];

opts = {'EdgeColor','none',...
        'FaceLighting','gouraud',...
        'AmbientStrength',0.5,...
        'DiffuseStrength',0.5,...
        'SpecularStrength',0.3,...
        'SpecularExponent',50,...
        'BackFaceLighting','reverselit'};
V = vdm.vals;
nv = length(V);
if isempty(vdm.clim)
    clim = prctile(V,[5,95]);
else
    clim = vdm.clim;
end
if vdm.logdisp
    V = real(log(V));
end
if nv == 0
    opts = [{'Facecolor','flat'},opts];
elseif nv == size(vdm.faces,1)
    opts = [{'FaceVertexCData',V,'Facecolor','flat'},opts];
elseif nv == size(vdm.vertices,1)
    opts = [{'FaceVertexCData',V,'Facecolor','interp'},opts];
else
    error('Invalid CData input length.');
end
opts = [{'Vertices',vdm.vertices,'Faces',vdm.faces},opts];

if nargin==3
    [fdir,fname,~] = fileparts(fname);
else
    fname = '';
end

% Generate figure, but don't save
if any(flag==0)
    h = genFig(1,clim,opts,vdm.label,vdm.logdisp);
end

% Save figure as JPEG:
if any(flag==1)
    if isempty(fname)
        [fname,fdir] = uiputfile('*.jpg','Save figure as:',fullfile(vdm.elxdir,'VDM.jpg'));
        if ~ischar(fname)
            fname = '';
        end
    else
        fname = [fname,'.jpg'];
    end
    if ~isempty(fname)
        h = genFig(2,clim,opts,vdm.label,vdm.logdisp);
        saveas(h.hfig,fullfile(fdir,fname));
    end
end

% Save object rotation movie:
if any(flag==2)
    if isempty(fname)
        [fname,fdir] = uiputfile('*.mp4','Save video as:',fullfile(vdm.elxdir,'VDM.mp4'));
        if ~ischar(fname)
            fname = '';
        end
    else
        fname = [fname,'.mp4'];
    end
    if ~isempty(fname)
        hh = genFig(1,clim,opts,vdm.label,vdm.logdisp);
        a = 0:4:359;
        v = VideoWriter(fullfile(fdir,fname),'MPEG-4');
        open(v);
        for i = 1:length(a)
            view(hh.plot(1).axes,90+a(i),0);
            lightangle(hh.plot(1).lightAngle,60+a(i),-30);
            writeVideo(v,getframe(hh.hfig));
        end
        close(v);
        h = [h,hh];
    end
end

% Save patch object as WRL:
if any(flag==3)
    if isempty(fname)
        [fname,fdir] = uiputfile('*.wrl','Save surface as:',fullfile(vdm.elxdir,'VDM.wrl'));
        if ~ischar(fname)
            fname = '';
        end
    else
        fname = [fname,'.wrl'];
    end
    if ~isempty(fname)
        if isempty(h)
            h = genFig(1,clim,opts,vdm.label,vdm.logdisp);
        end
        vrml(h.plot(1).axes,fullfile(fdir,fname),'noedgelines');
    end
end

% Save as X3D (XHTML):
if any(flag==4)
    if isempty(fname)
        [fname,fdir] = uiputfile('*.x3d','Save VDM as:',fullfile(vdm.elxdir,'VDM.x3d'));
        if ~ischar(fname)
            fname = '';
        end
    else
        fname = [fname,'.x3d'];
    end
    if ~isempty(fname)
        h = genFig(1,clim,opts,vdm.label,vdm.logdisp);
        disp('check');
        ha = h.plot.axes;
        data.tags.numobjects = 0;
        data.tags.xhtml = false;
        data.tags.folder = fdir;
        data.tags.filename = fname;
        data.tags.options = struct(...
            'output','both',...
            'width',600,'height',800,...
            'headlight',true,...
            'title','test',...
            'interactive',true);
        [data,loc_scene] = X3Dheader(data,[],ha);
        data = figurex3d(haxis,data,loc_scene);
        writex3dfile(folder,filename,data);
    end
end

function h = genFig(n,clim,opts,label,logdisp)
h.hfig = figure('Colormap',jet(128),'Position',[500 500 600*n 800],'Units','normalized');
for i = 1:n
    h.plot(i).axes = axes(h.hfig,'Position',[(i-1)/2 0 1/n 1],'CLim',clim);
    h.plot(i).lightAngle = lightangle(60+180*(i-1),-30);
    h.plot(i).patch = patch(h.plot(i).axes,opts{:});
    axis(h.plot(i).axes,'equal','off','tight');
    view(h.plot(i).axes,90+180*(i-1),0);
    if logdisp
        caxis(h.plot(i).axes,log(clim));
    end
end
h.cbar = colorbar(h.plot(1).axes,'FontSize',20,'AxisLocation','in',...
    'Position',[ 1/n + (n-2)*0.05 , 0.2 , 1/(n*40) , 0.6 ]);
yt = linspace(clim(1),clim(2),6);
if logdisp
%     yt = 0:0.2:3;
%     yclim = exp(h.cbar.YTick([1,end]));
%     yt(yt<yclim(1) | yt>yclim(2)) = [];
    h.cbar.YTick = log(yt);
    h.cbar.YTickLabel = yt;
else
    h.cbar.YTick = yt;
end
ht = title(h.cbar,['\fontsize{30}{0}\selectfont',label],'Interpreter','latex','Units','normalized','Position',[0.5,1.2,0]);
if n==1
    ht.HorizontalAlignment = 'right';
end


