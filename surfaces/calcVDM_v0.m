function h = calcVDM_v0(varargin)
% h = calcVDM(elxdir,dt,clim)
% h = calcVDM(mask,voxsz,parname,dt,clim)

switch nargin
    case 3
        fdir = varargin{1};
        dt = varargin{2};
        clim = varargin{3};
        [pout,pin] = readPtsFile(fdir);
        fv = load(fullfile(fdir,'fv.mat'));
    case 5
        
        parname = varargin{3};
        dt = varargin{4};
        clim = varargin{5};
        fdir = fileparts(parname);
        
        % Generate surface mesh:
        fv = mask2surf(varargin{1:2});
        
        % Save points:
        ptsname = fullfile(fdir,'inputPoints.txt');
        fid = fopen(ptsname,'w');
        nv = size(fv.vertices,1);
        fprintf(fid,'point\n%u\n',nv);
        for i = 1:nv
            fprintf(fid,'%f %f %f\n',fv.vertices(i,:));
        end
        fclose(fid);

        % Transform mesh points:
        str = strcat('/opt/X11/bin/xterm -geometry 170x50 -T "(Transformix)" -e ''',...
                     '/usr/local/bin/transformix -out "',fdir,'" -tp "',parname,...
                     '" -def "',ptsname,'"''');
        system(str);

        % Load transformed points:
        pin = fv.vertices;
        pout = readPtsFile(fdir);
        
    otherwise
        error('Invalid inputs.')
end

% Calculate relative change in area:
fv.method = 'dA';
fv.vertices = pout;
nf = size(fv.faces,1);
dA = zeros(nf,1);
fprintf('Calculating face areas ...\n');
for i = 1:nf
    A0 = area3D(pin(fv.faces(i,1),:),pin(fv.faces(i,2),:),pin(fv.faces(i,3),:));
    A1 = area3D(pout(fv.faces(i,1),:),pout(fv.faces(i,2),:),pout(fv.faces(i,3),:));
    dA(i) = A1/A0;
end
dA = dA.^(1/dt);

% Show map on most recent time point:
fv.facevertexcdata = dA;
save(fullfile(fdir,'fv.mat'),'-struct','fv');
h = surfmap(2,fv,clim);

function [pout,pin] = readPtsFile(elxdir)
disp('Reading points from file ...')
fid = fopen(fullfile(elxdir,'outputpoints.txt'),'r');
str = fread(fid,'*char')';
fclose(fid);
stro = regexp(str,'OutputPoint = \[ (\S+) (\S+) (\S+) ]','tokens');
stri = regexp(str,'InputPoint = \[ (\S+) (\S+) (\S+) ]','tokens');
nv = length(stro);
pout = zeros(nv,3);
pin = zeros(nv,3);
for i = 1:nv
    pout(i,:) = str2double(stro{i});
    pin(i,:) = str2double(stri{i});
end

function A = area3D(p1,p2,p3)
a = norm(p1-p2);
b = norm(p2-p3);
c = norm(p3-p1);
s = (a+b+c)/2;
A = sqrt(s*(s-a)*(s-b)*(s-c));