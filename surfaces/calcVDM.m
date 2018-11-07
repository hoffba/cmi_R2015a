function [vdm,h] = calcVDM_v2(segfname,elxdir,dt)
% [vdm,h] = calcVDM_v1(segfname,elxdir,dt)

h = [];
vdm = [];

if nargin==0
    [fname,fpath] = uigetfile('*.mhd','Select Segmentation:');
    if ischar(fname)
        segfname = fullfile(fpath,fname);
    else
        return;
    end
end
if nargin<2
    elxdir = uigetdir(pwd,'Select Elastix directory:');
    if ~ischar(elxdir)
        return;
    end
end
if nargin<3
    dt = str2double(inputdlg('Scale factor:','VDM Scale',1,{'1'}));
    if isempty(dt)
        return;
    elseif isnan(dt)
        error('Invalid scale factor.');
    end
end

units = 'yr';
vdm = struct(...
    'faces',{[]},'vertices',{[]},... % For plotting VDM
    'vertices_orig',{[]},...         % Vertices at baseline
    'map',{struct('method',{'dJ/dt','dSA/dt','GC','MC'},...
                  'label',{sprintf('$$\\dot{J}(/%s)$$',units),...
                           '$${\dot{A}(\%/yr)}$$',...
                           '$${GC}$$',...
                           '$${MC}$$'},...
                  'clim',{exp([-1,1]),[-100,100],[-0.1,0.1],[-0.3,0.3]},...
                  'logdisp',{true,false,false,false},...
                  'vals',{[],[],[],[]})},...
    'elxdir',{elxdir},...
    'scale',{dt},...
    'units',{units});

vdmfname = fullfile(elxdir,'vdm.mat');

% Find transform parameter file:
parname = dir(fullfile(elxdir,'TransformParameters.*.txt'));
ind = max(cellfun(@(x)str2double(x(21)),{parname(:).name}));
parname = sprintf('TransformParameters.%u.txt',ind);

% Find Jacobian rate:
JRfname = fullfile(elxdir,'spatialJacobianRate.mhd');
if exist(JRfname,'file')
    
else
end


Dname = fullfile(elxdir,'deformationField.mhd');
if ~exist(Dname,'file')
    fprintf('Generating deformation fields...\n');
    str = sprintf('cd %s ; /usr/local/bin/transformix -out ./ -tp ./%s -def all',elxdir,parname);
    [stat,cmdout] = system(str,'-echo');
    if stat
        error(cmdout);
    end
end

% Generate surface mesh:
if exist(segfname,'file')
    [mask,~,fov,info] = readMHD(segfname);
    d = size(mask);
    voxsz = fov./d;
else
    error('Could not find segmenation file: %s\n',segfname);
end
% Matrix space to real space:
vdm.direction = info.SliceOrient;
vdm.origin = info.SlicePos;
moff = (d-1).*voxsz/2;
n = cross(vdm.direction(1:3),vdm.direction(4:6));
M = [ reshape(vdm.direction,3,2),n',round(vdm.origin+moff,3)' ; zeros(1,3),1];

fprintf('Generating and transforming surface mesh ...\n');
[fv,x,y,z] = mask2surf(mask,voxsz);
vdm.vertices_orig = fv.vertices;
vdm.faces = fv.faces;
nv = size(fv.vertices,1);
V = M * [fv.vertices,ones(nv,1)]';
clear fv;
V = transformMesh(V(1:3,:)',fullfile(elxdir,parname));
V = (M \ [V';ones(1,nv)])';
vdm.vertices = V(:,1:3);

% Calculate VDM surface values:
fprintf('Reading deformation fields ...\n');
dX = readMHD(Dname);
fprintf('Calculating spatial Jacobian ...\n');
F = surfVal(def2jac(dX,voxsz),x,y,z,vdm.vertices_orig);
fprintf('Calculating spatial rate Jacobian ...\n');
FF = surfVal(def2jac(dX/dt,voxsz),x,y,z,vdm.vertices_orig);
np = size(F,1);
vdm.map(1).vals = zeros(np,1);
f = zeros(3); % J
ff = zeros(3);% dJ/dt
gi = mod(1:np,round(np/20))==0;
fprintf('Calculating Jacobian determinant rate ...\n');
for i = 1:np
    f(:) = F(i,:);
    ff(:) = FF(i,:);
    vdm.map(1).vals(i) = det(f) * trace( f \ ff )/3;
    if gi(i)
        fprintf('%u%% Complete\n',sum(gi(1:i))*5);
    end
end

% Save results:
save(vdmfname,'-struct','vdm');

% Display VDM and save JPEG:
[~,edirname] = fileparts(elxdir);
edirname = edirname(8:end);
h = savesurf(vdm,1,1,fullfile(elxdir,sprintf('%s_VDM_J.jpg',edirname)));

function V = surfVal(M,x,y,z,p)
nv = size(M,4);
V = zeros(size(p,1),nv);
for i = 1:nv
    V(:,i) = interp3(y,x,z,M(:,:,:,i),p(:,2),p(:,1),p(:,3));
end
