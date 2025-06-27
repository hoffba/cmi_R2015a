function fn_tfi = inverseTransform(fn_tf,fn_ref,fn_seg,info_hom)
% Generate inverse transform from registration
% Inputs:
%   fn_tf = file name of TransformParameters*.txt file defining forward reg
%   fn_ref = file name of registration reference image (e.g. EXP)
%   fn_seg = file name of segmentation associated with ref image
%   info_hom = Nifti info of homologous space (e.g. INS)

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
                     'SP_a',1,...
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
                          'ImageSampler',imsamp,...
                          'WriteResultImage','true',...
                          'SP_a',elx.Schedule{2}.SP_a/100);
end

% Make temp directory
elxdir = fileparts(fn_tf);
tdir = fullfile(elxdir,'temp_inv_elastix');
if ~isfolder(tdir)
    mkdir(tdir);
end

% Fix parameter files - make sure file paths are correct
fixTransformParameter(fn_tf);

% Generate Elastix call and run
cmdstr = elx.sysCmd(tdir,'wait',true,'f',fn_ref,'m',fn_ref,segin{:});
system(cmdstr);

% Prepare gometry variables
d = info_hom.ImageSize([2 1 3]);
voxsz = info_hom.PixelDimensions([2,1,3]);
orient = info_hom.Transform.T * diag([-1 -1 1 1]);
orient = orient([2,1,3,4],[2,1,3,4])'/diag([voxsz,1]);
orig = orient(1:3,4);
orient = reshape(orient(1:3,1:3),1,[]);

% Copy and rename TransformParameter files
fn_save = fullfile(elxdir,'InverseTransformParameters.0.txt');
txt = fileread(fullfile(tdir,'TransformParameters.0.txt'));
txt = regexprep(txt,'\(InitialTransformParametersFileName [^\)]*',...
    '(InitialTransformParametersFileName "NoInitialTransform"');
writelines(txt,fn_save);

fn_tfi = fn_save;
txt = fileread(fullfile(tdir,'TransformParameters.1.txt'));
txt = regexprep(txt,'\(InitialTransformParametersFileName [^\)]*',...
    ['\(InitialTransformParametersFileName "',regexprep(fn_tfi,'\\','\\\\'),'"']);
txt = regexprep(txt,'\(Size [^\)]*',...
    sprintf('(Size %f %f %f',d));
txt = regexprep(txt,'\(Spacing [^\)]*',...
    sprintf('(Spacing %f %f %f',voxsz));
txt = regexprep(txt,'\(Origin [^\)]*',...
    sprintf('(Origin %f %f %f',orig));
txt = regexprep(txt,'\(Direction [^\)]*',...
    sprintf('(Direction %f %f %f %f %f %f %f %f %f',orient));
writelines(txt,fullfile(elxdir,'InverseTransformParameters.1.txt'));

% Cleanup
fn = dir(fullfile(tdir,'*'));
for i = 3:numel(fn)
    delete(fullfile(tdir,fn(i).name));
end
rmdir(tdir);


