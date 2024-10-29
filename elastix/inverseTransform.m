function fn_tfi = inverseTransform(fn_tf,fn_ref,fn_seg)
% Generate inverse transform from registration
% Inputs:
%   fn_tf = file name of TransformParameters*.txt file defining forward reg
%   fn_ref = file name of registration reference image (e.g. EXP)
%   fn_seg = file name of segmentation associated with ref image

if nargin==3
    imsamp = 'RandomSparseMask';
    segin = {'fmask',fn_seg};
else
    imsamp = 'RandomCoordinate';
    segin = {};
end

% Find Elastix parameter file
elxdir = fileparts(fn_tf);
fn_elxpar = dir(fullfile(elxdir,'ElastixParameters_*.txt'));
elxchk = isempty(fn_elxpar);
if ~elxchk
    fn_elxpar = fullfile(elxdir,fn_elxpar(end).name);
end

% Initialize elastix inputs
elx = ElxClass;
elx.setTx0(fn_tf);
elx.addStep('affine','Metric','DisplacementMagnitudePenalty',...
                     'ImageSampler',imsamp,...
                     'NumberOfResolutions',1,...
                     'NumberOfSpatialSamples',5000,...
                     'MaximumNumberOfIterations',1000,...
                     'SP_a',10,...
                     'SP_A',50,...
                     'FixedImagePyramidSchedule',[1 1 1],...
                     'MovingImagePyramidSchedule',[1 1 1]);
if elxchk
    elx.addStep('warp','Metric','DisplacementMagnitudePenalty',...
                       'ImageSampler',imsamp)
else
    elx.loadPar(fn_elxpar);
    elx.setPar(2,'Metric',{'DisplacementMagnitudePenalty'},...
                          'Metric0Weight',1,...
                          'GridSpacingSchedule',[5 5 5 1 1 1],...
                          'ImageSampler',imsamp,...
                          'WriteResultImage','true',...
                          'SP_a',elx.Schedule{2}.SP_a/10,...
                          'NumberOfSpatialSamples',[5000,10000],...
                          'MaximumNumberOfIterations',[10000,5000]);
end

% Make temp directory
elxdir = fileparts(fn_tf);
tdir = fullfile(elxdir,'temp_inv_elastix');
if ~isfolder(tdir)
    mkdir(tdir);
end

% Generate Elastix call and run
cmdstr = elx.sysCmd(tdir,'wait',true,'f',fn_ref,'m',fn_ref,segin{:});
system(cmdstr);

% Copy and rename TransformParameter files
fn_save = fullfile(elxdir,'InverseTransformParameters.0.txt');
txt = fileread(fullfile(tdir,'TransformParameters.0.txt'));
txt = regexprep(txt,'\(InitialTransformParametersFileName [^\)]*',...
    '(InitialTransformParametersFileName "NoInitialTransform"');
writelines(txt,fn_save);

fn_tfi = regexprep(fn_save,'\\','\\\\');
txt = fileread(fullfile(tdir,'TransformParameters.1.txt'));
txt = regexprep(txt,'\(InitialTransformParametersFileName [^\)]*',...
    ['\(InitialTransformParametersFileName "',fn_tfi,'"']);
writelines(txt,fullfile(elxdir,'InverseTransformParameters.1.txt'));

% Cleanup
fn = dir(fullfile(tdir,'*'));
for i = 3:numel(fn)
    delete(fullfile(tdir,fn(i).name));
end
rmdir(tdir);


