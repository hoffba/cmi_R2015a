function info = init_niftiinfo(label,voxsz,dtype,d)

info = struct( 'Description',label,...
               'PixelDimensions',voxsz,...
               'Datatype',dtype,...
               'ImageSize',d(1:3),...
               'Version','NIfTI1',...
               'Qfactor',-1,...
               'SpaceUnits','Millimeter',...
               'TimeUnits','None',...
               'SliceCode','Unknown',...
               'FrequencyDimension',0,...
               'PhaseDimension',0,...
               'SpatialDimension',3);