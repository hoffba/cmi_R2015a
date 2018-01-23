function patch2wrl(h)


[fname,fdir] = uiputfile('*.wrl','Save 3D Surface:');
if ischar(fname)
    vrml(ha1,fullfile(fdir,fname),'noedgelines');
end
