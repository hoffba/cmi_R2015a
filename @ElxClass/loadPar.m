% ElxClass function
% Load parameters from file
function stat = loadPar(self,fname)

stat = false;
if (nargin==1) || isempty(fname)
    [fname,path] = uigetfile('*.txt','Load ElastixParameters.txt',...
                             'MultiSelect','on');
    if iscell(fname) || ischar(fname)
        fname = strcat(path,fname);
    else
        fname = {};
    end
elseif ~(ischar(fname) || iscellstr(fname))
    fname = {};
end
if ischar(fname)
    fname = {fname};
end

% Loop over files to load:
nf = length(fname);
newpars = cell(1,nf);
for i = 1:nf
    if exist(fname{i},'file')
        s = readElxTxt2Struct(fname{i});
        if all(isfield(s,{'Registration',... % Look for all required Elastix Pars
                          'Metric',...
                          'ImageSampler',...
                          'Interpolator',...
                          'ResampleInterpolator',...
                          'Resampler',...
                          'Transform',...
                          'Optimizer',...
                          'FixedImagePyramid',...
                          'MovingImagePyramid'}))
            newpars{i} = s;
            stat = true;
        else
            error('Invalid parameter file.')
        end
        % Set single-value inputs to Nres vectors for RegClass
        n = newpars{i}.NumberOfResolutions;
        if length(newpars{i}.SP_A)==1
            newpars{i}.SP_A = newpars{i}.SP_A * ones(1,n);
        end
        if length(newpars{i}.SP_a)==1
            newpars{i}.SP_a = newpars{i}.SP_a * ones(1,n);
        end
        if length(newpars{i}.NumberOfSpatialSamples)==1
            newpars{i}.NumberOfSpatialSamples = ...
                newpars{i}.NumberOfSpatialSamples * ones(1,n);
        end
        if length(newpars{i}.MaximumNumberOfIterations)==1
            newpars{i}.MaximumNumberOfIterations = ...
                newpars{i}.MaximumNumberOfIterations * ones(1,n);
        end
    end
end
newpars(cellfun(@isempty,newpars)) = [];
if ~isempty(newpars)
    self.Schedule = [self.Schedule,newpars];
    stat = true;
end