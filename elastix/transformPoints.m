function tpts = transformPoints(fn_tf,pts)
% Inputs:
%   fn_tf = file name for transform
%   pts   = point you'd like to transform, [x,y,z]
%       ** pts need to be in spatial coordinates determined by image
%          orientation and voxel size

procdir = fileparts(fn_tf);
np = size(pts);

% Write points to file
fn_pts = fullfile(procdir,'inputPoints.txt');
fid = fopen(fn_pts,'w');
fprintf(fid,'point\n%d\n',size(pts,1));
fprintf(fid,'%f %f %f\n',pts');
fclose(fid);

% Transformix command
cmdstr = sprintf('transformix -def "%s" -out "%s" -tp "%s"',fn_pts,procdir,fn_tf);
system(cmdstr);

% Return transformed points
tpts = nan(np(1),3);
fid = fopen(fullfile(procdir,'outputpoints.txt'),'r');
for i = 1:np(1)
    str = fgetl(fid);
    str = extractBetween(str,'OutputPoint = [ ',' ]');
    tpts(i,:) = double(string(strsplit(str{1})));
end