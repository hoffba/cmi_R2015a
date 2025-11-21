function mapYACTA_AirwaysResults(ID,ydir,fn_airways)
% Maps YACTA airway results to airway surface using nearest value
% Inputs:
%       ydir = path to yacta processing directory, containing AirwayResults.mhd
%       fn_airways = path to airways file

% Read airways
fprintf('Reading airways: %s\n',fn_airways);
[A,label,fov,orient,info,fnameOut] = cmi_load(1,[],fn_airways);
d = size(A);

% Read YACTA results
fn_results = fullfile(ydir,'AirwayResults.mhd');
fprintf('Reading YACTA results: %s\n',fn_results)
R = readMHD(fn_results);
[dR(1),dR(2),dR(3),nv] = size(R);

% Check for validity of input dimensions
if ~all(d==dR)
    error('Dimensions of airways and results do not match.');
end

% Generate airways surface triangulation
fprintf('Triangulating airway surface\n')
voxsz = fov./d;
[X,Y,Z] = meshgrid( (1:d(2)) * voxsz(2),...
                    (1:d(1)) * voxsz(1),...
                    (1:d(3)) * voxsz(3));
[faces, vertices] = isosurface(X,Y,Z,logical(A),0.5);
TR = triangulation(faces, vertices);
clear faces vertices

% find coordinates to resuts skeleton
skel = any(R,4);
sk_ind = find(skel);
[sky,skx,skz] = ind2sub(d,sk_ind);
sk_coord = [skx,sky,skz] .* voxsz([2,1,3]);

% Plot airways
hf = figure;
ht = trisurf(TR,'EdgeColor','none');
ha = ht.Parent;
axis(ha,'equal');
set(ha,'XColor','none','YColor','none','ZColor','none',...
       'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{},...
       'XGrid','off','YGrid','off','ZGrid','off');
view(ha,0,0);

% Find face colors
fprintf('Finding face colors\n')
nF = length(ht.Faces);
C = nan(nF,nv);
for i = 1:nF
    G = mean(ht.Vertices(ht.Faces(i,:),:),1);
    sk_d = sqrt(sum((G-sk_coord).^2,2));
    [~,ind] = min(sk_d);
    C(i,:) = squeeze(R(sky(ind),skx(ind),skz(ind),:))';
end

% Loop over results vectors
if nv==4
    label = ["Lumen Area","Wall Area","Wall Percentage","Bronchiectasis Index"];
else
    label = string(1:nv);
end
map = jet(128);
colormap(hf,map);
colorbar
for i = 1:nv
    nC = C(:,i);
    mxC = max(nC);
    if ~mxC
        mxC = 1;
    end
    nC = nC/max(nC)*128;
    clim(ha,[0,mxC]);
    ht.FaceVertexCData = squeeze(ind2rgb(uint8(nC),map));
    title(label(i));
    saveas(hf,fullfile(ydir,strcat(ID,'_',matlab.lang.makeValidName(label(i)),'.tif')));
end

delete(hf);
