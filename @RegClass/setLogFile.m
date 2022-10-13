function setLogFile(self,fname)

if ischar(fname) && endsWith(fname,'.txt')
    self.fn_log = fname;
end