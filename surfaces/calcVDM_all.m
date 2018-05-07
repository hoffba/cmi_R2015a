function [vdm,h] = calcVDM(varargin)
% [vdm,h] = calcVDM(elxdir,'Name',Value,...)
% [vdm,h] = calcVDM(vdm,elxdir,'Name',Value,...)
%
% Methods:
%   1 = linearly scaled J : V = (J-1)/N
%   2 = exponentially scaled J : J^(1/N)
%   3 = Surface tangent J^(1/N)
%   4 = dJ/dN
%   5 = dA
%
% Inputs:
%   vdm = VDM results structure
%   elxdir = elxreg_ directory
%   'Name'/Value options [default]:
%       'scale'     [1]
%       'method'    [1]
%       'mask'
%       'voxsz'     [ones(1,3)]
%       'clim'      [auto]
%       'logdisp'	[true]
%       'units'     ['yr']
%       'faces'
%       'vertices'
% Outputs:
%   vdm = struct with fields:
%       .mask       = binary image segmentation
%       .voxsz      = voxel dimensions
%       .faces      = surface mesh faces
%       .vertices   = surface mesh vertices
%       .vals       = surface color values
%       .units      = string describing normalization units (yr, mmHg, etc.)
%       .elxdir     = elxreg_ directory
%       .scale      = normalization scale factor
%       .method     = type of VDM map to generate
%       .clim       = color limits for display
%       .logdisp    = T/F display colors in log scale?

[vdm,flag,x,y,z] = parseinputs(varargin{:});
nv = length(vdm.vertices);

% Generate surface map based on method:
switch vdm.method
    case 1 % linear |J|
        vdm.label = sprintf('$$({J}-1)/{%s}+1$$',vdm.units);
        vdm.vals = (surfVal(loadData(vdm.elxdir,'jdet'),x,y,z,vdm.vertices)-1)/vdm.scale + 1;
    case 2 % exponential |J|
        vdm.label = sprintf('$${J}^{(1/%s)}$$',vdm.units);
        vdm.vals = surfVal(loadData(elxdir,'jdet'),x,y,z,vdm.vertices).^(1/vdm.scale);
    case 3 % surface tangent |J|
        vdm.label = sprintf('$$Tangent({J})^{(1/%s)}$$',vdm.units);
        M = permute(reshape(surfVal(loadData(elxdir,'jmat'),x,y,z,vdm.vertices)',3,3,[]),[2,1,3]);

        % Determine planar Jacobian from isonormals:
        n = isonormals(x,y,z,mask,fv.vertices);
        vdm.vals = zeros(nv,1);
        for i = 1:size(n,1)
            R = vrrotvec2mat(vrrotvec([0,1,0],n(i,:)));
            tjac = R * M(:,:,i) * R';
            vdm.vals(i) = det(tjac(1:2,1:2));
        end
        clear M n;
        vdm.vals = vdm.vals.^(1/vdm.scale);

    case 4 % dJ/dt
        vdm.label = sprintf('$$\\dot{J}(%s^{-1})$$',vdm.units);
        
        dX = loadData(elxdir,'def');
        fprintf('Calculating spatial Jacobian ...\n');
        F = surfVal(def2jac(dX,voxsz),x,y,z,vdm.vertices(:,:,1));
        fprintf('Calculating spatial rate Jacobian ...\n');
        FF = surfVal(def2jac(dX/N,voxsz),x,y,z,vdm.vertices(:,:,1));

        np = size(F,1);
        vdm.vals = zeros(np,1);
        f = zeros(3); % J
        ff = zeros(3);% dJ/dt
        gi = mod(1:np,round(np/20))==0;
        fprintf('Calculating Jacobian rate determinant ...\n');
        for i = 1:np
            f(:) = F(i,:);
            ff(:) = FF(i,:);
            vdm.vals(i) = det(f) * trace( ff / f )/3;
            if gi(i)
                fprintf('%u%% Complete\n',sum(gi(1:i))*5);
            end
        end
        
    case 5 % dA
        
        vdm.label = sprintf('$${\\delta{Area}(\\%%^{%s^{-1}})}$$',vdm.units);
        
        % Save points:
        ptsname = fullfile(elxdir,'inputPoints.txt');
        fid = fopen(ptsname,'w');
        nv = size(fv.vertices,1);
        fprintf(fid,'point\n%u\n',nv);
        for i = 1:nv
            fprintf(fid,'%f %f %f\n',fv.vertices(i,:));
        end
        fclose(fid);
        
        p = loadData(elxdir,'pts');
        vdm.vertices = p;
        nf = size(vdm.faces,1);
        vdm.vals = zeros(nf,1);
        fprintf('Calculating face areas ...\n');
        for i = 1:nf
            vdm.vals(i) = area3D(p(vdm.faces(i,:),:,:));
        end
        vdm.vals = (vdm.vals-1)/vdm.scale;
end
clear x y z;

% Save resulting maps:
svname = fullfile(vdm.elxdir,'VDM.mat');
[fname,fdir] = uiputfile('*.mat','Save VDM as:',svname);
if ischar(fname)
    save(fullfile(fdir,fname),'-struct','vdm');
end

% Generate VDM figure:
h = savesurf(vdm,flag);

function [vdm,flag,x,y,z] = parseinputs(varargin)
x=[];y=[];z=[];
if isstruct(varargin{1}) && all(isfield(varargin{1},{'mask','voxsz','vertices','faces'}))
    vdm = varargin{1};
    varargin(1) = [];
else
    vdm = struct('mask',[],'voxsz',[],'vertices',[],'faces',[]);
end
p = inputParser;
addRequired(p,'elxdir',@(x)ischar(x)&&exist(x,'dir'));
addParameter(p,'scale',1,@isscalar);
addParameter(p,'method',1,@isscalar);
addParameter(p,'mask',[],@(x)islogical(x));
addParameter(p,'voxsz',[],@(x)isvector(x)&&(numel(x)==3));
addParameter(p,'clim',[],@(x)isvector(x)&&(numel(x)==2));
addParameter(p,'logdisp',true,@islogical);
addParameter(p,'units','yr',@ischar)
addParameter(p,'faces',[],@islogical);
addParameter(p,'vertices',[],@islogical);
addParameter(p,'disp',1,@(x)all(ismember(x,0:4)));
parse(p,varargin{:});
if isempty(vdm.vertices) && ismember('mask',p.UsingDefaults) && any(ismember({'faces','vertices'},p.UsingDefaults))
    error('Cannot calculate VDM without segmentation surface.');
end
if ~isempty(p.Results.mask)
    if isempty(vdm.vertices)
        answer = 'Yes';
    else
        answer = questdlg('Replace existing surface mesh?','Yes','No','Yes');
    end
    if strcmp(answer,'Yes')
        vdm.mask = p.Results.mask;
        if ~ismember('voxsz',p.UsingDefaults)
            vdm.voxsz = p.Results.voxsz;
        elseif isempty(vdm.voxsz)
            warning('Voxel dimensions not input, using ones.');
            vdm.voxsz = ones(1,3);
        end
        [fv,x,y,z] = mask2surf(vdm.mask,vdm.voxsz);
        vdm.faces = fv.faces;
        vdm.vertices = fv.vertices;
    end
end
vdm.elxdir = p.Results.elxdir;
vdm.scale = p.Results.scale;
vdm.method = p.Results.method;
vdm.clim = p.Results.clim;
vdm.logdisp = p.Results.logdisp;
vdm.units = p.Results.units;
flag = p.Results.disp;


function V = surfVal(M,x,y,z,p)
nv = size(M,4);
V = zeros(size(p,1),nv);
for i = 1:nv
    V(:,i) = interp3(x,y,z,M(:,:,:,i),p(:,1),p(:,2),p(:,3));
end

function A = loadData(elxdir,opt)
C = {'def'  ,'deformationField.mhd'     ,'-def all';...
     'jdet' ,'SpatialJacobian.mhd'      ,'-jac all';...
     'jmat' ,'FullSpatialJacobian.mhd'  ,'-jacmat all';...
     'pts'  ,'outputpoints.txt'         ,'-def "./inputPoints.txt"'};
i = find(strcmp(opt,C(:,1)),1);
if isempty(i)
    error('Invalid input.')
end
fname = fullfile(elxdir,C{i,2});
if ~exist(fname,'file')
    % Need to generate from transform:
    parname = dir(fullfile(elxdir,'TransformParameters.*.txt'));
    ind = max(cellfun(@(x)str2double(x(21)),{parname(:).name}));
    str = sprintf(['/opt/X11/bin/xterm -geometry 170x50 -T "(Transformix)" -e ''',...
             'cd %s ; /usr/local/bin/transformix -out "./" -tp ',...
             '"./TransformParameters.%u.txt" %s'''],elxdir,ind,C{i,3});
    system(str);
end
if i==4
    A = readPtsFile(fname);
else
    A = readMHD(fname);
end

function p = readPtsFile(fname)
% Vertex points: [nv x (x/y/z) x (in/out)]
disp('Reading points from file ...')
fid = fopen(fname,'r');
str = fread(fid,'*char')';
fclose(fid);
pat = ' = \[ (\S+) (\S+) (\S+) ]';
tok = cellfun(@(x)str2double(x)',regexp(str,['InputPoint',pat],'tokens'),'UniformOutput',false);
p = [tok{:}]';
tok = cellfun(@(x)str2double(x)',regexp(str,['OutputPoint',pat],'tokens'),'UniformOutput',false);
p(:,:,2) = [tok{:}]';

function A = area3D(p)
A = [0,0];
for i = 1:2
    a = norm(p(1,:,i)-p(2,:,i));
    b = norm(p(2,:,i)-p(3,:,i));
    c = norm(p(3,:,i)-p(1,:,i));
    s = (a+b+c)/2;
    A(i) = sqrt(s*(s-a)*(s-b)*(s-c));
end
A = A(2)/A(1);
    