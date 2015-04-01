function stat = isDICOM(fname)
% Detect whether a file is in DICOM format
% Output: fid = true - it is a DICOM
%               false - it is not a DICOM
%               NaN - could not open file

stat = NaN;
if ischar(fname)
    fname = {fname};
end
if nargin && iscellstr(fname)
    nf = numel(fname);
    stat = NaN(1,nf);
    for i = 1:nf
        fid = fopen(fname{i},'r');
        if fid>0
            fseek(fid,128,'bof');
            tag = fread(fid,4,'char=>char')';
            stat(i) = strcmp(tag,'DICM')';
            if ~stat(i)
                fseek(fid,0,'bof');
                tag = fread(fid,2,'uint32');
                stat(i) = (isequal(tag,[8 0]) || isequal(tag,[134217728 0]) || ...
                        isequal(tag,[8 4]) || isequal(tag,[134217728 67108864]));
            end
            fclose(fid);
        end
    end
end