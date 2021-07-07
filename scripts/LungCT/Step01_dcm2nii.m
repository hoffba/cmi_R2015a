function data = Step01_dcm2nii(procdir, fnames, tags)

% Use dicom catalog csv file to extract information.
% home_pwd = home directory
% ID_list = current ID being processed
% list = files associated with ID (dicom Exp and dicom Ins)

if numel(fnames) ~= numel(tags)
    error('Inputs fnames and tags must have same number of elements')
end
% outFormat = '.nii.gz';
data = struct('tag',tags,'img',{[],[]},'voi',{[],[]});
nf = length(fnames);
for i=1:nf
    if ~isfolder(procdir)
        mkdir(procdir)
    end
    
    [img,~,~,orient,info] = readDICOM(fnames{i});
    data(i).img.mat = img;
    
    % Generate NIfTI info:
    data(i).img.info = struct('Filename',fullfile(procdir,sprintf('%s_%s_%u.nii',...
                                    info.meta.PatientID,info.meta.StudyDate,info.meta.SeriesNumber)),...
                              'Description',info.meta.StudyDescription,...
                              'PixelDimensions',[info.meta.PixelSpacing',info.meta.SlcThk],...
                              'Datatype','int16',...
                              'ImageSize',size(img),...
                              'Version','NIfTI1',...
                              'Qfactor',-1,...
                              'SpaceUnits','Millimeter',...
                              'TimeUnits','None',...
                              'SliceCode','Unknown',...
                              'FrequencyDimension',0,...
                              'PhaseDimension',0,...
                              'SpatialDimension',3,...
                              'TransformName','Sform',...
                              'Transform',affine3d(orient'));
    
%% Commented by BAH 2021-06-28 in favor of using other DICOM reader
%     dicm2nii(fnames{i}, procdir, outFormat);
%     list_file = dir(fullfile(procdir,'*.gz'));
    
%     if ~isempty(list_file)
%         fprintf('Reading image from file (#%u): %s\n',i,list_file.name);
%         str_file = fullfile(list_file.folder,list_file.name);
% 
%         data(i).img.info = niftiinfo(str_file);
%         data(i).img.mat = niftiread(data(i).img.info);
%         
%         data(i).voi.info = data(i).img.info;
%         data(i).voi.info.Datatype = 'int8';
%         data(1).voi.info.BitsPerPixel = 8;
%         
%         delete(str_file)
%     end
end