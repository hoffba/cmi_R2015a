function check_EI = CTlung_ReadImages(regObj,fnames,dcmnames)

imchk  = cellfun(@(x)exist(x,'file'),fnames);
if ~all(imchk)
    check_EI = true;
end

fprintf('\nReading image data from file ... ID = %s\n',ID)
for i = 1:2
    if imchk(i)
        fprintf('   ... from NiFTi\n');
        loadname = fnames{i};
    else
        fprintf('   ... from DICOM\n');
        loadname = dcmnames{i};
    end
    regObj.cmiObj(i).loadImg(0,loadname);
end

