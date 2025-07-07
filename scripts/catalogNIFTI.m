function T = catalogNIFTI(searchdir)

T = [];

try

    fn = dir(fullfile(searchdir,'**','*.nii.gz'));
    
    vars = {...
            'Tag',               'cellstr';...
            'UMlabel',           'cellstr';...
            'SeriesDescription', 'cellstr';...
            'PatientName',       'cellstr';...
            'StudyDate',         'cellstr';...
            'SeriesNumber',      'single';...
            'ConvolutionKernel', 'cellstr';...
            'SliceThickness',    'single';...
            'Slices',            'single';...
            'DataType',          'cellstr';...
            'DataPath',          'cellstr'};
    nv = size(vars,1);
    Tdef = table('Size',[1,nv],'VariableTypes',vars(:,2)','VariableNames',vars(:,1)');
    T = table('Size',[0,nv],'VariableTypes',vars(:,2)','VariableNames',vars(:,1)');
    for i = 1:numel(fn)
        t = Tdef;
        fname = fn(i).name;
        t.UMlabel = {extractBefore(fname,'_')};
        t.PatientName = {extractBefore(fname,'_')};
        t.DataType = {'NIFTI'};
        t.DataPath = {fullfile(fn(i).folder,fname)};
        tag = '';
        if contains(fname,'exp','IgnoreCase',true)
            tag = 'Exp';
        end
        if contains(fname,'ins','IgnoreCase',true)
            tag = 'Ins';
        end
        if contains(fname,'segmentation','IgnoreCase',true)
            tag = [tag,'Label'];
        end
        if contains(fname,'warped','IgnoreCase',true)
            tag = 'InsReg';
        end
        t.Tag = {tag};
        T = [T;t];
    end

catch err
    disp(getReport(err))
end