function [ success ] = writeCSV( filename, A, outlabels)
%UNTITLED Summary of this function goes here
% Inputs:
%   A = matrix to write out
%   outlabels = set of cells containing text for labels
%
%   Detailed explanation goes here


% Open file
fid = fopen(filename,'w');
if fid==-1
    fprintf('****ERROR opening %s\n',filename);
    success = 0;
    return;
end

% write labels
for ii=1:length(outlabels)
    fprintf(fid,'%s,',outlabels{ii});
end; 
fprintf(fid,'\n');

% write array
for ii=1:size(A,1)
    for jj=1:size(A,2)
        fprintf(fid,'%d,',A(ii,jj));
    end
    fprintf(fid,'\n');
end
fclose(fid);

success = 1;

end

