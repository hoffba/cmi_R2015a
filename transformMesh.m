function V = transformMesh(V,tpname)

[tdir,tname] = fileparts(tpname);
if isempty(tdir)
    error('Input transform parameter file with full directory path.');
end

nv = size(V,1);

% Save points:
fprintf('Saving mesh points ...\n');
ptsname = fullfile(tdir,'inputPoints.txt');
fid = fopen(ptsname,'w');
fprintf(fid,'point\n%u\n',nv);
for i = 1:nv
    fprintf(fid,'%f %f %f\n',V(i,:));
end
fclose(fid);

% Generate call to transformix:
str = sprintf('cd %s ; /usr/local/bin/transformix -out ./ -tp ./%s.txt -def inputPoints.txt',tdir,tname);
[stat,cmdout] = system(str,'-echo');
if stat
    error(cmdout);
end

% Load transformed points:
opstr = fullfile(tdir,'outputpoints.txt');
if exist(opstr,'file')
    fprintf('Loading transformed points ...\n');
    fid = fopen(opstr,'r');
    str = fread(fid,'*char')';
    fclose(fid);
    stro = regexp(str,'OutputPoint = \[ (\S+) (\S+) (\S+) ]','tokens');
    nv = length(stro);
    for i = 1:nv
        V(i,:) = str2double(stro{i});
    end
else
    error('File outputpoints.txt not found.');
end
