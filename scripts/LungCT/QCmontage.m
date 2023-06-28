function cdata = QCmontage(tag,im_bg,im_over,voxsz,fname)
% Generate figures for displaying montage overlays with segmentations or PRM

nf = size(im_over,3);

% Find system screen size
set(0,'units','pixels');
scnsz = get(0,'screensize');

% Initialize figure:
hf = figure('visible','off','Position',scnsz,'PaperPositionMode','auto');
ha = axes(hf,'DataAspectRatio',voxsz);
im = imagesc(im_bg(:,:,1),[-1000 0]);
axis image
set(ha,'Units','pixels');
npos = plotboxpos(ha);
npos(1:2) = npos(1:2)+1;
npos(3:4) = npos(3:4)-2;
F = getframe(hf,npos);
montsz = [round(sqrt(nf)),nan];

% Generate frames:
switch tag
    case {'seg','reg'}
        colormap(hf,'gray');
        hold(ha,'on');
        hroi = plot(ha,1,1,'*m','MarkerSize',10);
        hold(ha,'off');
        for i = 1:nf
            im.CData = im_bg(:,:,i);
            [voir,voic] = find(edge(im_over(:,:,i),'Canny'));
            set(hroi,'XData',voic,'YData',voir);
            F(i) = getframe(hf,npos);
        end
    case 'prm'
        colormap(hf,[ gray(256) ;...
                      0 1 0 ;... % Norm
                      1 1 0 ;... % fSAD
                      1 0 0 ;... % Emph
                      1 0 1 ;... % PD
                      1 1 1 ]);  % NS
        mat = (max(min(im_bg,0),-1000)+1000)*256/1000 - 0.5;
        ind = im_over>0;
        mat(ind) = im_over(ind) + 255.5;
        ha.CLim = [0 261];
        for i = 1:nf
            im.CData = mat(:,:,i);
            F(i) = getframe(hf,npos);
        end
end

% Generate montage
mat = zeros(size(F(1).cdata,1),size(F(1).cdata,2),3,nf);
for i = 1:nf
    [im,map] = frame2im(F(i));
    if ~isempty(map)
        im = ind2rgb(im,map);
    end
    mat(:,:,:,i) = double(im);
end
figure(hf)
montage(mat/255,'Size',montsz,'Parent',ha);
[~,tname] = fileparts(fname);
t = title(tname,'Interpreter','none');

% Fit figure to axes
axpos = tightPosition(ha);
ha.Position(1:2) = ha.Position(1:2)-axpos(1:2);
hf.Position(3:4) = [axpos(3:4)] - 1 + [0 t.Extent(4)];

% Print the figure:
pipeline_save_fig(hf,fname);

cdata = print(hf,'-RGBImage','-noui');

delete(hf);
